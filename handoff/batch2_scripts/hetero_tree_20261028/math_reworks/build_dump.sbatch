#!/bin/bash
#SBATCH --job-name=mrw_build
#SBATCH --partition=norm
#SBATCH --constraint=x6342
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
# math reworks (branch claude/tree-math-reworks off claude/tree-dist 9005a8c77): build S into V, then the fast gate
set -uo pipefail
: ${S:?} ${V:?}
G=/mnt/nasapps/production/gcc/15.2.0; PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=$G/bin:/mnt/nasapps/production/cmake/4.0.3/bin:/mnt/nasapps/production/python/3.12/bin:$PATH
export LD_LIBRARY_PATH=$G/lib64:${LD_LIBRARY_PATH:-}
mkdir -p $V; cd $V; date > build_started.txt; hostname > build_host.txt
cmake -S $S -B $V/build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DUSE_OPENMP_OFFLOAD=OFF -DUSE_ARCHOPT=ON -DCMAKE_INSTALL_PREFIX=$V/install -DCMAKE_C_COMPILER=$G/bin/gcc -DCMAKE_CXX_COMPILER=$G/bin/g++ -DCMAKE_Fortran_COMPILER=$G/bin/gfortran -DPython3_EXECUTABLE=$PY -DBLA_VENDOR=OpenBLAS > cmake.log 2>&1 || { echo cmake failed; exit 6; }
cmake --build $V/build --parallel 32 > make.log 2>&1 || { grep -n -i error -B2 -A6 make.log | head -60; exit 7; }
cmake --install $V/build > install.log 2>&1 && date > build_finished.txt
bash $S/scripts/run_fast_gate.sh $V/build > fast_gate.log 2>&1; echo $? > fast_gate.exit
echo "build done, fast gate exit $(cat fast_gate.exit)"; grep -E "tests passed|tests failed" fast_gate.log
