#!/bin/bash
# merge_sm: the merge branch tip plus the prob_state_marg patch, built apart from merge_master2 (its install serves running jobs)
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_sm
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:/mnt/nasapps/production/python/3.12/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd $B/build
hostname > $B/build_host.txt; date > $B/build_started.txt
cmake $B/src -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$B/install -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_ARCHOPT=ON -DUSE_OPENMP=ON > $B/cmake.log 2>&1 && make -j 32 install > $B/make.log 2>&1
echo $? > $B/build.exit; date > $B/build_finished.txt
[ "$(cat $B/build.exit)" = 0 ] || exit 1
cat > $B/fast_gate.sbatch <<EOS
#!/bin/bash
#SBATCH -J fastgate_sm
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=$B/fast_gate.out
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:/mnt/nasapps/production/python/3.12/bin:\$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd $B/src && scripts/run_fast_gate.sh $B/build > $B/fast_gate.log 2>&1; echo \$? > $B/fast_gate.exit
EOS
sbatch --parsable $B/fast_gate.sbatch > $B/fast_gate.jobid
