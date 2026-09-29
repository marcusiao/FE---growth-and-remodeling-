"""Read-only comparison of completed matched short cycles; writes a new report."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
try:
    import h5py
except ImportError:
    # Read-only fallback to the HDF5 libraries already installed with DOLFIN.
    import ctypes as ct
    import ctypes.util
    class NativeFile:
        def __init__(self,path,mode):
            assert mode=='r'
            self.lib=ct.CDLL(ctypes.util.find_library('hdf5_openmpi'))
            self.hl=ct.CDLL(ctypes.util.find_library('hdf5_openmpi_hl'))
            self.lib.H5Fopen.argtypes=[ct.c_char_p,ct.c_uint,ct.c_longlong];self.lib.H5Fopen.restype=ct.c_longlong
            self.lib.H5Fclose.argtypes=[ct.c_longlong]
            self.hl.H5LTget_dataset_ndims.argtypes=[ct.c_longlong,ct.c_char_p,ct.POINTER(ct.c_int)]
            self.hl.H5LTget_dataset_info.argtypes=[ct.c_longlong,ct.c_char_p,ct.POINTER(ct.c_ulonglong),ct.POINTER(ct.c_int),ct.POINTER(ct.c_size_t)]
            self.hl.H5LTread_dataset_double.argtypes=[ct.c_longlong,ct.c_char_p,ct.POINTER(ct.c_double)]
            self.handle=self.lib.H5Fopen(str(path).encode(),0,0)
            assert self.handle>=0
        def __enter__(self):return self
        def __exit__(self,*args):self.lib.H5Fclose(self.handle)
        def __getitem__(self,key):
            key=key.encode();n=ct.c_int()
            assert self.hl.H5LTget_dataset_ndims(self.handle,key,ct.byref(n))>=0
            dims=(ct.c_ulonglong*n.value)();kind=ct.c_int();size=ct.c_size_t()
            assert self.hl.H5LTget_dataset_info(self.handle,key,dims,ct.byref(kind),ct.byref(size))>=0
            a=np.empty(tuple(dims),dtype=np.float64)
            assert self.hl.H5LTread_dataset_double(self.handle,key,a.ctypes.data_as(ct.POINTER(ct.c_double)))>=0
            if kind.value==0:
                assert np.all(a==np.floor(a)) and np.max(abs(a))<2**53
                a=a.astype(np.int64)
            return a
    class h5py:
        File=NativeFile

def fields(folder,phase):
    root=ET.parse(folder/(phase+'.xdmf')).getroot()
    grids=root.findall('./Domain/Grid/Grid');out={}
    with h5py.File(folder/(phase+'.h5'),'r') as h:
        def read(item):
            filename,key=item.text.strip().split(':',1)
            assert filename==phase+'.h5'
            return h[key][...]
        xyz=read(grids[0].find('./Geometry/DataItem'))
        cells=read(grids[0].find('./Topology/DataItem'))
        for g in grids:
            t=float(g.find('Time').attrib['Value'])
            out[t]={a.attrib['Name']:read(a.find('DataItem')) for a in g.findall('Attribute')}
    return xyz,cells,out

def main():
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path)
    args=p.parse_args();a=json.loads((args.baseline/'summary.json').read_text());b=json.loads((args.candidate/'summary.json').read_text())
    assert a['status']==b['status']=='COMPLETED_EXPLORATORY_NOT_VALIDATED'
    assert len(a['states'])==len(b['states'])==24
    for key in ('evolution','materials','IOP_MPa','CSFP_MPa','mesh_manifest','runtime_mesh_manifest'):
        assert a[key]==b[key],key
    assert a['settings']==b['settings']
    report=dict(meaning='Settings sensitivity relative to q4/Newton, not absolute accuracy',
                baseline=str(args.baseline),candidate=str(args.candidate),phases={},state_differences=[],fields={})
    report['total_seconds']={'baseline':a['elapsed_seconds'],'candidate':b['elapsed_seconds'],
                             'speedup':a['elapsed_seconds']/b['elapsed_seconds']}
    for phase in ('loading','remodeling','unloading'):
        sa=[s for s in a['states'] if s['phase']==phase];sb=[s for s in b['states'] if s['phase']==phase]
        report['phases'][phase]={}
        for label,states in (('baseline',sa),('candidate',sb)):
            report['phases'][phase][label]=dict(steps=len(states),mean_s=float(np.mean([s['solver']['seconds'] for s in states])),
                 mean_iterations=float(np.mean([s['solver']['iterations'] for s in states])))
        report['phases'][phase]['candidate']['timing_sums']={k:sum(s['solver']['timing'][k] for s in sb) for k in sb[0]['solver']['timing']}
        report['phases'][phase]['candidate']['counts']={k:sum(s['solver']['counts'][k] for s in sb) for k in sb[0]['solver']['counts']}
        xa,ca,fa=fields(args.baseline,phase);xb,cb,fb=fields(args.candidate,phase)
        assert np.array_equal(xa,xb) and np.array_equal(ca,cb) and fa.keys()==fb.keys()
        report['fields'][phase]={}
        for t in fa:
            v={};tags=fa[t]['Material_ID'].ravel();assert np.array_equal(tags,fb[t]['Material_ID'].ravel())
            for name in ('Displacement','Theta_Jg','J_total','E1_total','von_Mises'):
                v[name]={}
                for region,tag in (('whole',None),('PT',9),('LC',8)):
                    mask=np.arange(len(fa[t][name])) if tag is None else (np.unique(ca[tags==tag]) if name=='Displacement' else np.flatnonzero(tags==tag))
                    va=fa[t][name][mask];vb=fb[t][name][mask];diff=vb-va;den=np.linalg.norm(va)
                    v[name][region]=dict(max_abs=float(np.max(abs(diff))),rms=float(np.sqrt(np.mean(diff**2))),
                         relative_l2=float(np.linalg.norm(diff)/den) if den>1e-14 else None)
            report['fields'][phase][str(t)]=v
    for sa,sb in zip(a['states'],b['states']):
        assert (sa['phase'],sa['step'],sa['load_factor'])==(sb['phase'],sb['step'],sb['load_factor'])
        assert sb['sampled_J_min']>0 and sb['sampled_Je_min']>0
        if sb['solver']:assert sb['solver']['residual']<=sb['solver']['tolerance']
        row=dict(phase=sa['phase'],step=sa['step'])
        for k in ('PT_volume','LC_volume','PT_mean_theta','theta_min'):
            row[k]=dict(baseline=sa[k],candidate=sb[k],difference=sb[k]-sa[k],percent_difference=100*(sb[k]/sa[k]-1))
        report['state_differences'].append(row)
    dest=args.candidate/'comparison_to_q4.json'
    with dest.open('x') as f:json.dump(report,f,indent=2,allow_nan=False)
    lines=['# q3/BFGS versus q4/Newton short-cycle comparison','',report['meaning'],'',
           f"Total elapsed: baseline {a['elapsed_seconds']:.2f}s; candidate {b['elapsed_seconds']:.2f}s; speedup {a['elapsed_seconds']/b['elapsed_seconds']:.3f}x.",
           '','| Phase | Baseline mean solve (s) | Candidate mean solve (s) | Solve speedup |',
           '|---|---:|---:|---:|']
    for phase,v in report['phases'].items():
        old=v['baseline']['mean_s'];new=v['candidate']['mean_s']
        lines.append(f'| {phase} | {old:.3f} | {new:.3f} | {old/new:.3f} |')
    lines+=['','## Phase endpoint differences','',
            'Displacement relative L2 is an unweighted nodal vector norm; volume differences are relative to the baseline at the same step. Near-zero reference fields have no relative norm.','',
            '| Phase | PT volume difference (%) | LC volume difference (%) | PT displacement relative L2 (%) | LC displacement relative L2 (%) |',
            '|---|---:|---:|---:|---:|']
    for phase in report['phases']:
        row=[x for x in report['state_differences'] if x['phase']==phase][-1]
        end=max(report['fields'][phase],key=float);f=report['fields'][phase][end]['Displacement']
        lines.append(f"| {phase} | {row['PT_volume']['percent_difference']:.6g} | {row['LC_volume']['percent_difference']:.6g} | {100*f['PT']['relative_l2']:.6g} | {100*f['LC']['relative_l2']:.6g} |")
    lines+=['','Both settings changed together; this test cannot isolate their individual effects. Timings are historical rather than simultaneous controlled benchmarks. BFGS action timing includes base linear solving; do not sum both. Neither run establishes absolute accuracy, locking resolution or full-trajectory robustness. Full simulation remains on hold.']
    with (args.candidate/'COMPARISON.md').open('x') as f:f.write('\n'.join(lines)+'\n')
    print(dest);print(json.dumps(report['total_seconds']));print(json.dumps(report['phases'],indent=2))

if __name__=='__main__':main()

