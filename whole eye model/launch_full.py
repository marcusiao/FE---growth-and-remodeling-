"""User-selected strict q3/BFGS full cycle, detached with persistent log."""
from pathlib import Path
from datetime import datetime
import hashlib,json,os,subprocess,sys

root=Path(__file__).resolve().parent
passed=root/'results/20260924_143754_479346'
r=json.loads((passed/'summary.json').read_text())
assert r['status']=='COMPLETED_EXPLORATORY_NOT_VALIDATED' and len(r['states'])==24
assert r['settings']['smoke'] and r['settings']['e1_threshold']==.005
assert r['settings']['pressure_mode']=='reference' and r['setup']['quadrature_degree']==3
for s in r['states']:
    assert s['sampled_J_min']>0 and s['sampled_Je_min']>0
    if s['solver']:assert s['solver']['residual']<=s['solver']['tolerance']
for name in ('model_whole_eye_hex.py','bfgs_solver.py','hex_geometry.py'):
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==hashlib.sha256((passed/'source'/name).read_bytes()).hexdigest(),name
for p in Path('/proc').iterdir():
    if not p.name.isdigit():continue
    try:args=(p/'cmdline').read_bytes().split(b'\0')
    except OSError:continue
    if any(Path(a.decode(errors='replace')).name=='model_whole_eye_hex.py' for a in args):
        raise RuntimeError('Another whole-eye model is running: '+p.name)
folder=root/'launches'/datetime.now().strftime('%Y%m%d_%H%M%S_%f');folder.mkdir(parents=True)
command=[sys.executable,'-u',str(root/'model_whole_eye_hex.py'),'--full-cycle','--e1-threshold','0.005','--pressure-mode','reference']
with (folder/'simulation.log').open('wb') as log:
    proc=subprocess.Popen(command,cwd=root,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
data=dict(status='LAUNCHED_NOT_COMPLETED',pid=proc.pid,command=command,passed_smoke=str(passed),
          log=str(folder/'simulation.log'),schedule=[100,500,100],candidate='q3_BFGS_strict_36point_checks')
(folder/'launch.json').write_text(json.dumps(data,indent=2));print(json.dumps(data,indent=2))
