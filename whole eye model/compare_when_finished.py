"""Finish only the authorized comparison; never launches a simulation."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path);p.add_argument('pid',type=int)
a=p.parse_args();status=a.candidate/'comparison_status.json'
def save(state,**extra):
    tmp=status.with_suffix('.tmp');tmp.write_text(json.dumps(dict(status=state,**extra),indent=2));tmp.replace(status)
save('WAITING_FOR_SHORT_TEST')
try:
    while True:
        r=json.loads((a.candidate/'summary.json').read_text())
        if r['status']=='COMPLETED_EXPLORATORY_NOT_VALIDATED':break
        if r['status']=='FAILED':raise RuntimeError('Short test failed; inspect its saved traceback')
        if not Path('/proc',str(a.pid)).exists():raise RuntimeError('Short test ended without completion')
        time.sleep(20)
    result=subprocess.run([sys.executable,str(Path(__file__).with_name('compare_runs.py')),str(a.baseline),str(a.candidate)],capture_output=True,text=True)
    (a.candidate/'comparison.log').write_text(result.stdout+'\n'+result.stderr)
    if result.returncode:raise RuntimeError('Comparison failed; see comparison.log')
    save('COMPARISON_COMPLETED',report=str(a.candidate/'comparison_to_q4.json'))
except Exception as exc:
    save('FAILED',error=str(exc));raise
