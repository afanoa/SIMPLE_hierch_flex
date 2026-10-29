#!/usr/bin/env python3
"""Deliver raw FLEX state maps with an existing SIMPLE executable in one allocation.

The input contract pins every scientific input. Each task owns a copied project;
only a complete, validated map set is published. No project is merged or edited
outside task scratch. Execution requires Slurm; --plan only validates the contract.
"""
import argparse
import concurrent.futures
import datetime
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import signal
import socket
import struct
import subprocess
import sys
import threading
import time


SUFFIXES = ('', '_even_unfil', '_odd_unfil')
TERMINAL = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY',
            'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE'}


def sha256(path, offset=0):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        stream.seek(offset)
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def read_contract(path):
    contract = json.loads(Path(path).read_text())
    required = {'version', 'executable', 'inputs', 'project', 'stage_files',
                'nstates', 'mskdiam', 'pgrp', 'box'}
    if not required <= contract.keys() or contract['version'] != 1:
        raise ValueError('Expected version 1 input contract with all required fields')
    for key in ('nstates', 'box'):
        if type(contract[key]) is not int or contract[key] < 1:
            raise ValueError(key + ' must be a positive integer')
    if not 0 < float(contract['mskdiam']) < float('inf'):
        raise ValueError('mskdiam must be finite and positive')
    if not re.fullmatch(r'[cCdD][1-9][0-9]*|[tToOiI]', contract['pgrp']):
        raise ValueError('Unsupported point-group spelling')
    inputs = contract['inputs']
    if not isinstance(inputs, dict) or not inputs:
        raise ValueError('inputs must map absolute paths to SHA256 digests')
    for name, digest in inputs.items():
        if not Path(name).is_absolute() or not re.fullmatch(r'[0-9a-f]{64}', digest):
            raise ValueError('Invalid pinned input: ' + name)
    staged = contract['stage_files']
    if not isinstance(staged, list) or not staged or len(set(staged)) != len(staged):
        raise ValueError('stage_files must be a unique, nonempty list')
    if len({Path(p).name for p in staged}) != len(staged):
        raise ValueError('Stage file basenames collide')
    if contract['project'] not in staged or not contract['project'].endswith('.simple'):
        raise ValueError('project must be a staged .simple file')
    if contract['executable'] not in inputs or any(p not in inputs for p in staged):
        raise ValueError('Executable and staged files must be pinned inputs')
    if 'state_weights.txt' not in {Path(p).name for p in staged}:
        raise ValueError('The frozen FLEX state weight manifest must be staged')
    expected = contract.get('expected_payloads', {})
    names = {map_name(s, suffix) for s in range(1, contract['nstates'] + 1) for suffix in SUFFIXES}
    if expected and (set(expected) != names or any(not re.fullmatch(r'[0-9a-f]{64}', v) for v in expected.values())):
        raise ValueError('expected_payloads must pin every full/even_unfil/odd_unfil map')
    return contract


def map_name(state, suffix):
    return 'recvol_state%02d%s.mrc' % (state, suffix)


def make_plan(contract, parallel, workers, threads, cpus):
    if parallel not in ('yes', 'no') or min(workers, threads, cpus) < 1:
        raise ValueError('Invalid concurrency settings')
    active = min(workers, contract['nstates']) if parallel == 'yes' else 1
    required = active * threads + 1
    if required > cpus:
        raise ValueError(f'Need {required} CPUs: {active} workers x {threads} threads + 1 controller; have {cpus}')
    if parallel == 'no' and workers != 1:
        raise ValueError('--workers must be 1 when --parallel=no')
    plan = {'contract': contract, 'parallel': parallel, 'workers': active,
            'launcher_sha256': sha256(__file__),
            'threads_per_worker': threads, 'controller_cpus': 1, 'reserved_cpus': required,
            'allocation_cpus': cpus, 'descendant_scheduler_jobs': 0,
            'program': 'simple_exec prg=reconstruct3D m_estimator=flex rec_backend=gridding',
            'publication': 'raw full/even_unfil/odd_unfil maps only; source project unchanged',
            'postprocess': 'no' if parallel == 'yes' else 'yes'}
    encoded = json.dumps(plan, sort_keys=True, separators=(',', ':')).encode()
    plan['fingerprint'] = hashlib.sha256(encoded).hexdigest()
    return plan


def verify_inputs(contract, workers):
    def verify(item):
        path, expected = item
        if sha256(path) != expected:
            raise ValueError('Input checksum changed: ' + path)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(verify, contract['inputs'].items()))


def mrc_record(path, box):
    with Path(path).open('rb') as stream:
        header = stream.read(1024)
    if len(header) != 1024 or header[208:212] != b'MAP ':
        raise ValueError('Invalid MRC header: ' + str(path))
    integers = struct.unpack('<256i', header)
    floats = struct.unpack('<256f', header)
    if integers[:4] != (box, box, box, 2) or integers[23] < 0:
        raise ValueError('MRC grid/mode mismatch: ' + str(path))
    offset = 1024 + integers[23]
    if Path(path).stat().st_size != offset + box**3 * 4:
        raise ValueError('Incomplete MRC payload: ' + str(path))
    frame = list(integers[:10]) + list(floats[10:16]) + list(integers[16:19]) + list(floats[49:52])
    if not all(math.isfinite(v) for v in frame):
        raise ValueError('Nonfinite MRC frame: ' + str(path))
    return {'file': Path(path).name, 'bytes': Path(path).stat().st_size,
            'sha256': sha256(path), 'payload_sha256': sha256(path, offset), 'frame': frame}


def validate_maps(directory, states, contract):
    result = []
    for state in states:
        for suffix in SUFFIXES:
            record = mrc_record(Path(directory) / map_name(state, suffix), contract['box'])
            expected = contract.get('expected_payloads', {}).get(record['file'])
            if expected is not None and record['payload_sha256'] != expected:
                raise ValueError('Payload differs from the reference: ' + record['file'])
            result.append(record)
    if any(r['frame'] != result[0]['frame'] for r in result):
        raise ValueError('State/half map frames differ')
    return result


class Runner:
    def __init__(self, plan, output, job_id):
        self.plan, self.output, self.job_id = plan, Path(output), job_id
        self.stop = threading.Event()
        self.mutex = threading.Lock()
        self.processes = {}
        self.launch_pending = 0
        self.worker_local = threading.local()
        self.cpu_slots = queue.SimpleQueue()
        self.taskset = None
        if hasattr(os, 'sched_getaffinity'):
            available = sorted(os.sched_getaffinity(0))
            if len(available) < plan['reserved_cpus']:
                raise ValueError('Worker groups exceed available CPU affinity')
            self.taskset = shutil.which('taskset')
            if self.taskset is None:
                raise RuntimeError('Linux execution requires taskset for disjoint worker CPU groups')
            width = plan['threads_per_worker']
            for worker in range(plan['workers']):
                self.cpu_slots.put(available[1+worker*width:1+(worker+1)*width])
        self.environment = {k: v for k, v in os.environ.items() if not k.startswith('SIMPLE_')}
        build = Path(plan['contract']['executable']).parent.parent
        self.environment.update(SIMPLE_PATH=str(build), SIMPLE_QSYS='local', SIMPLE_EMAIL='nil', SIMPLE_SEED='1',
                                OMP_NUM_THREADS=str(plan['threads_per_worker']), OPENBLAS_NUM_THREADS='1',
                                MKL_NUM_THREADS='1', OMP_THREAD_LIMIT=str(plan['threads_per_worker']))
        self.environment.update(OMP_PLACES='cores', OMP_PROC_BIND='close')
        self.environment['PATH'] = str(build / 'bin') + os.pathsep + self.environment.get('PATH', '')

    def assign_worker_cpus(self):
        self.worker_local.cpus = self.cpu_slots.get() if self.taskset else []

    def record_processes(self):
        atomic_json(self.output / 'processes.json', {'job_id': self.job_id, 'host': socket.gethostname(),
                    'pids': sorted(self.processes), 'launch_pending': self.launch_pending, 'updated_utc': utc_now()})

    def run_command(self, command, directory, label):
        if self.stop.is_set():
            raise RuntimeError('Cancelled before launch')
        if self.taskset:
            command = [self.taskset, '--cpu-list', ','.join(map(str,self.worker_local.cpus)), *command]
        atomic_json(directory / (label + '_command.json'), command)
        start = time.monotonic()
        with (directory / (label + '.log')).open('wb') as log:
            with self.mutex:
                self.launch_pending += 1
                self.record_processes()
            try:
                process = subprocess.Popen(command, cwd=directory, env=self.environment,
                                           stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            except BaseException:
                with self.mutex:
                    self.launch_pending -= 1
                    self.record_processes()
                raise
            try:
                with self.mutex:
                    self.processes[process.pid] = process
                    self.launch_pending -= 1
                    self.record_processes()
                while process.poll() is None:
                    if self.stop.wait(0.1):
                        raise RuntimeError('Worker cancelled')
                if process.returncode:
                    raise RuntimeError(f'{label} failed with exit {process.returncode}; see {directory}')
            finally:
                if process.poll() is None:
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                with self.mutex:
                    self.processes.pop(process.pid, None)
                    self.record_processes()
        return time.monotonic() - start

    def task(self, states):
        name = 'all_states' if self.plan['parallel'] == 'no' else 'state_%02d' % states[0]
        task_root = self.output / name
        task_root.mkdir(exist_ok=True)
        contract = self.plan['contract']
        saved = task_root / 'complete.json'
        if saved.exists():
            result = json.loads(saved.read_text())
            if result['fingerprint'] != self.plan['fingerprint']:
                raise ValueError('Task checkpoint belongs to another plan')
            current = validate_maps(task_root / result['attempt'], states, contract)
            if current != result['maps']:
                raise ValueError('Completed task output changed: ' + name)
            return result
        attempt = 1
        while (task_root / ('attempt_%03d' % attempt)).exists():
            attempt += 1
        directory = task_root / ('attempt_%03d' % attempt)
        directory.mkdir()
        start = time.monotonic()
        for source in contract['stage_files']:
            if self.stop.is_set():
                raise RuntimeError('Cancelled while staging')
            destination = directory / Path(source).name
            shutil.copyfile(source, destination)
            if sha256(destination) != contract['inputs'][source]:
                raise ValueError('Staged input checksum changed: ' + source)
        project = Path(contract['project']).name
        exe = contract['executable']
        update_s = self.run_command([exe, 'prg=update_project', 'projfile=' + project,
                                     'qsys_name=local', 'mkdir=no', 'nthr=1'], directory, 'update_project')
        command = [exe, 'prg=reconstruct3D', 'projfile=' + project, 'mkdir=no', 'oritype=ptcl3D',
                   'nstates=' + str(contract['nstates']), 'm_estimator=flex', 'rec_backend=gridding',
                   'mskdiam=' + str(contract['mskdiam']), 'pgrp=' + contract['pgrp'], 'objfun=euclid',
                   'qsys_name=local', 'nthr=' + str(self.plan['threads_per_worker']),
                   'postprocess=' + self.plan['postprocess']]
        if self.plan['parallel'] == 'yes':
            command.append('state=' + str(states[0]))
        native_s = self.run_command(command, directory, 'reconstruct')
        if 'SIMPLE_RECONSTRUCT3D NORMAL STOP' not in (directory / 'reconstruct.log').read_text(errors='replace'):
            raise ValueError('Native normal-stop marker missing: ' + name)
        result = {'fingerprint': self.plan['fingerprint'], 'task': name, 'attempt': directory.name,
                  'states': states, 'maps': validate_maps(directory, states, contract),
                  'worker_cpus': self.worker_local.cpus,
                  'update_s': update_s, 'native_s': native_s, 'task_wall_s': time.monotonic() - start}
        atomic_json(saved, result)
        return result

    def run(self):
        tasks = [[s] for s in range(1, self.plan['contract']['nstates'] + 1)]
        if self.plan['parallel'] == 'no':
            tasks = [list(range(1, self.plan['contract']['nstates'] + 1))]
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.plan['workers'], initializer=self.assign_worker_cpus) as pool:
            futures = [pool.submit(self.task, states) for states in tasks]
            try:
                results = [future.result() for future in concurrent.futures.as_completed(futures)]
            except BaseException:
                self.stop.set()
                for future in futures:
                    future.cancel()
                raise
        return sorted(results, key=lambda r: r['states'][0])


def reject_live_previous(output, current_job):
    path = output / 'processes.json'
    if not path.exists():
        return
    previous = json.loads(path.read_text())
    if not previous['pids'] and not previous.get('launch_pending', 0):
        return
    if previous['host'] == socket.gethostname() and not previous.get('launch_pending', 0):
        for pid in previous['pids']:
            try:
                os.killpg(pid, 0)
            except ProcessLookupError:
                continue
            raise RuntimeError('Prior worker PID may still be alive; refusing to overlap: ' + str(pid))
        return
    job = previous['job_id']
    if job == current_job or not re.fullmatch(r'[0-9]+(?:_[0-9]+)?', job):
        raise RuntimeError('Cannot prove prior workers are stopped on the other host')
    status = subprocess.check_output(['sacct', '-X', '-j', job, '-n', '-P', '-o', 'State'], text=True)
    states = [line.split('|')[0].split()[0].rstrip('+') for line in status.splitlines() if line.strip()]
    if not states or any(state not in TERMINAL for state in states):
        raise RuntimeError('Prior allocation is not confirmed terminal')


def execute(plan, output, resume, job_id):
    output = Path(output).absolute()
    output.mkdir(parents=True, exist_ok=resume)
    with (output / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan_file = output / 'plan.json'
        if plan_file.exists():
            if json.loads(plan_file.read_text()) != plan:
                raise ValueError('Resume plan/input/concurrency differs from the saved plan')
        else:
            atomic_json(plan_file, plan)
        reject_live_previous(output, job_id)
        start = time.monotonic()
        verify_inputs(plan['contract'], min(plan['allocation_cpus'] - 1, 8))
        published = output / 'maps'
        if published.exists():
            if not resume:
                raise FileExistsError(published)
            complete = json.loads((published / 'manifest.json').read_text())
            records = validate_maps(published, range(1, plan['contract']['nstates'] + 1), plan['contract'])
            if complete['fingerprint'] != plan['fingerprint'] or complete['maps'] != records:
                raise ValueError('Published output does not match its manifest')
            return complete
        runner = Runner(plan, output, job_id)
        old_handlers = {}
        if threading.current_thread() is threading.main_thread():
            def stop(signum, frame):
                runner.stop.set()
                raise InterruptedError('Controller received signal ' + str(signum))
            for signum in (signal.SIGTERM, signal.SIGINT):
                old_handlers[signum] = signal.signal(signum, stop)
        try:
            results = runner.run()
            verify_inputs(plan['contract'], min(plan['allocation_cpus'] - 1, 8))
            temporary = output / 'maps.pending'
            if temporary.exists():
                shutil.rmtree(temporary)
            temporary.mkdir()
            records = []
            for task in results:
                for record in task['maps']:
                    source = output / task['task'] / task['attempt'] / record['file']
                    destination = temporary / record['file']
                    shutil.copyfile(source, destination)
                    with destination.open('rb') as stream:
                        os.fsync(stream.fileno())
                    if sha256(destination) != record['sha256']:
                        raise ValueError('Output changed during publication')
                    records.append(record)
            if any(record['frame'] != records[0]['frame'] for record in records):
                raise ValueError('State/half map frames differ across workers')
            complete = {'status': 'COMPLETE', 'fingerprint': plan['fingerprint'], 'maps': records,
                        'tasks': results, 'reserved_cpus': plan['reserved_cpus'], 'job_id': job_id,
                        'wall_s_including_validation': time.monotonic() - start, 'finished_utc': utc_now(),
                        'sources_unchanged': True, 'reference_payloads_checked': bool(plan['contract'].get('expected_payloads')),
                        'scope': plan['publication']}
            atomic_json(temporary / 'manifest.json', complete)
            sync_directory(temporary)
            os.rename(temporary, published)
            sync_directory(output)
            atomic_json(output / 'status.json', {'status': 'COMPLETE', 'finished_utc': utc_now()})
            return complete
        except BaseException as error:
            runner.stop.set()
            atomic_json(output / 'status.json', {'status': 'FAILED', 'error': str(error), 'finished_utc': utc_now()})
            raise
        finally:
            for signum, handler in old_handlers.items():
                signal.signal(signum, handler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--parallel', choices=('no', 'yes'), default='no')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--cpus', type=int, required=True, help='Total allocation budget, including one controller CPU')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--plan', action='store_true', help='Print the validated plan; no files, hashes, or processes')
    args = parser.parse_args()
    plan = make_plan(read_contract(args.contract), args.parallel, args.workers, args.threads, args.cpus)
    if args.plan:
        print(json.dumps(plan, indent=2))
        return
    job = os.environ.get('SLURM_JOB_ID', '')
    allocated = int(os.environ.get('SLURM_CPUS_PER_TASK', '0'))
    if not re.fullmatch(r'[0-9]+(?:_[0-9]+)?', job) or args.cpus > allocated:
        parser.error('Execution requires one Slurm allocation with SLURM_CPUS_PER_TASK >= --cpus')
    if hasattr(os, 'sched_getaffinity') and args.cpus > len(os.sched_getaffinity(0)):
        parser.error('Requested budget exceeds this process CPU affinity')
    result = execute(plan, args.output, args.resume, job)
    print(json.dumps({k: result[k] for k in ('status', 'reserved_cpus', 'wall_s_including_validation',
                                          'reference_payloads_checked')}, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError) as error:
        print('state-maps: ' + str(error), file=sys.stderr)
        sys.exit(1)
