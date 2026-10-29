#!/usr/bin/env python3
"""Offline orchestration tests; no SIMPLE executable or scientific computation."""
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import hetero_tree_state_maps as subject


FAKE = '''#!/usr/bin/env python3
import json,os,pathlib,struct,sys,time
args=dict(arg.split('=',1) for arg in sys.argv[1:])
if args['prg']=='update_project':
 pathlib.Path(args['projfile']).write_text('private project updated')
 sys.exit(0)
states=[int(args['state'])] if 'state' in args else list(range(1,int(args['nstates'])+1))
if os.environ.get('FAKE_FAIL_STATE') in [str(s) for s in states]:
 time.sleep(.3)
 sys.exit(7)
time.sleep(.03)
for state in states:
 for suffix in ('','_even_unfil','_odd_unfil'):
  header=bytearray(1024)
  struct.pack_into('<4i',header,0,2,2,2,2)
  struct.pack_into('<3i',header,28,2,2,2)
  struct.pack_into('<6f',header,40,4.,4.,4.,90.,90.,90.)
  struct.pack_into('<3i',header,64,1,2,3)
  if os.environ.get('FAKE_FRAME_STATE')==str(state):
   struct.pack_into('<f',header,196,1.)
  header[208:212]=b'MAP '
  payload=struct.pack('<8f',*([float(state)]*8))
  pathlib.Path('recvol_state%02d%s.mrc'%(state,suffix)).write_bytes(header+payload)
print('SIMPLE_RECONSTRUCT3D NORMAL STOP')
'''


class StateMapWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        exe = self.root/'build'/'bin'/'simple_exec'
        exe.parent.mkdir(parents=True)
        exe.write_text(FAKE)
        exe.chmod(0o755)
        project = self.root/'f.simple'
        project.write_text('original project')
        weights = self.root/'state_weights.txt'
        weights.write_text('frozen weights')
        self.contract = {'version':1,'executable':str(exe),'project':str(project),
                         'stage_files':[str(project),str(weights)],
                         'inputs':{str(p):subject.sha256(p) for p in (exe,project,weights)},
                         'nstates':2,'mskdiam':4,'pgrp':'c1','box':2,
                         'expected_payloads':{subject.map_name(s,x):hashlib.sha256(struct.pack('<8f',*([float(s)]*8))).hexdigest()
                                              for s in (1,2) for x in subject.SUFFIXES}}

    def tearDown(self):
        self.temporary.cleanup()

    def plan(self, mode='yes'):
        return subject.make_plan(self.contract,mode,2 if mode=='yes' else 1,2,5)

    def test_budget_reserves_controller_and_default_serial(self):
        with self.assertRaisesRegex(ValueError,'Need 5 CPUs'):
            subject.make_plan(self.contract,'yes',2,2,4)
        plan=self.plan('no')
        self.assertEqual(plan['reserved_cpus'],3)
        self.assertEqual(plan['postprocess'],'yes')
        self.assertEqual(self.plan()['postprocess'],'no')
        with self.assertRaises(ValueError):
            subject.make_plan(self.contract,'no',2,2,8)

    def test_serial_parallel_payloads_and_private_projects(self):
        serial=subject.execute(self.plan('no'),self.root/'serial',False,'123')
        parallel=subject.execute(self.plan(),self.root/'parallel',False,'123')
        self.assertEqual(serial['maps'],parallel['maps'])
        self.assertEqual(Path(self.contract['project']).read_text(),'original project')
        self.assertEqual(len(parallel['maps']),6)
        commands=list((self.root/'parallel').glob('state_*/attempt_*/reconstruct_command.json'))
        self.assertEqual(len(commands),2)
        for p in commands:
            args=json.loads(p.read_text())
            self.assertIn('postprocess=no',args)
            self.assertFalse(any(a.startswith('nparts=') for a in args))
        self.assertFalse((self.root/'parallel'/'maps.pending').exists())

    def test_failure_keeps_publication_absent_and_resume_reuses_success(self):
        output=self.root/'failed'
        with patch.dict(os.environ,{'FAKE_FAIL_STATE':'2'}):
            with self.assertRaises(RuntimeError):
                subject.execute(self.plan(),output,False,'123')
        self.assertFalse((output/'maps').exists())
        self.assertEqual(json.loads((output/'processes.json').read_text())['pids'],[])
        result=subject.execute(self.plan(),output,True,'123')
        self.assertEqual(result['status'],'COMPLETE')
        self.assertEqual(result['tasks'][0]['attempt'],'attempt_001')
        self.assertEqual(result['tasks'][1]['attempt'],'attempt_002')

    def test_complete_resume_validates_outputs_and_pins(self):
        output=self.root/'resume'
        result=subject.execute(self.plan(),output,False,'123')
        self.assertEqual(result,subject.execute(self.plan(),output,True,'123'))
        changed=self.plan(); changed['threads_per_worker']=1
        with self.assertRaisesRegex(ValueError,'differs'):
            subject.execute(changed,output,True,'123')
        target=output/'maps'/subject.map_name(1,'')
        target.write_bytes(target.read_bytes()[:-4])
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            subject.execute(self.plan(),output,True,'123')

    def test_cross_worker_frame_mismatch_prevents_publication(self):
        output=self.root/'frame_mismatch'
        with patch.dict(os.environ,{'FAKE_FRAME_STATE':'2'}):
            with self.assertRaisesRegex(ValueError,'frames differ across workers'):
                subject.execute(self.plan(),output,False,'123')
        self.assertFalse((output/'maps').exists())

    def test_disjoint_cpu_groups_reserve_one_controller_cpu(self):
        with patch.object(subject.os,'sched_getaffinity',return_value={10,11,12,13,14},create=True), \
             patch.object(subject.shutil,'which',return_value='/usr/bin/taskset'):
            runner=subject.Runner(self.plan(),self.root,'123')
            runner.assign_worker_cpus(); first=runner.worker_local.cpus
            runner.assign_worker_cpus(); second=runner.worker_local.cpus
        self.assertEqual(first,[11,12])
        self.assertEqual(second,[13,14])

    def test_input_mutation_fails_before_processes(self):
        Path(self.contract['project']).write_text('changed source')
        with self.assertRaisesRegex(ValueError,'Input checksum changed'):
            subject.execute(self.plan(),self.root/'input_changed',False,'123')
        self.assertFalse((self.root/'input_changed'/'processes.json').exists())

    def test_contract_rejects_collisions_and_incomplete_reference(self):
        path=self.root/'contract.json'
        path.write_text(json.dumps(self.contract))
        self.assertEqual(subject.read_contract(path),self.contract)
        bad=dict(self.contract); bad['expected_payloads']={'missing.mrc':'0'*64}
        path.write_text(json.dumps(bad))
        with self.assertRaisesRegex(ValueError,'every'):
            subject.read_contract(path)

    def test_prior_live_process_group_and_unknown_launch_block_resume(self):
        output=self.root/'old'; output.mkdir()
        child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],start_new_session=True)
        try:
            subject.atomic_json(output/'processes.json',{'job_id':'123','host':socket.gethostname(),
                                                        'pids':[child.pid],'launch_pending':0})
            with self.assertRaisesRegex(RuntimeError,'still be alive'):
                subject.reject_live_previous(output,'123')
        finally:
            os.killpg(child.pid,signal.SIGTERM); child.wait()
        subject.reject_live_previous(output,'123')
        subject.atomic_json(output/'processes.json',{'job_id':'123','host':socket.gethostname(),
                                                    'pids':[],'launch_pending':1})
        with self.assertRaisesRegex(RuntimeError,'Cannot prove'):
            subject.reject_live_previous(output,'123')

    def test_cli_execution_requires_slurm_but_plan_does_not(self):
        path=self.root/'contract.json'; path.write_text(json.dumps(self.contract))
        command=[sys.executable,str(Path(subject.__file__)), '--contract',str(path),'--output',str(self.root/'cli'),
                 '--threads','2','--cpus','3']
        env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_')}
        planned=subprocess.run(command+['--plan'],env=env,capture_output=True,text=True)
        self.assertEqual(planned.returncode,0,planned.stderr)
        self.assertFalse((self.root/'cli').exists())
        run=subprocess.run(command,env=env,capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('requires one Slurm allocation',run.stderr)


if __name__=='__main__':
    unittest.main()
