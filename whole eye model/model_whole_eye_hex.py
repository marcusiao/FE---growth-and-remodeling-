"""Experimental Q1 whole-eye kinematic model with a verified legacy-DOLFIN mesh ordering."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import csv
import hashlib
import json
import shutil
import time
import traceback
from datetime import datetime
from pathlib import Path
import numpy as np
from hex_geometry import volumes, sampled_J

ROOT=Path(__file__).resolve().parent
MESH=ROOT.parents[2]/'01_mesh_conversion/whole_eye_hex_no_bm_20260923/mesh'
RUNTIME_MESH=ROOT.parent/'reordered_mesh'
# E in MPa. Tags4-7 are explicit exploratory placeholders, NOT literature values.
MATERIALS={1:('Sclera',3.),2:('Choroid',.2),4:('Retina',.1),5:('ON',.1),
           6:('Pia',.3),7:('Dura',3.),8:('LC',.3),9:('PT',.1)}
NU=.49

def save(path,value):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False));temp.replace(path)

def evolution(theta,e1,tags,crit):
    new=theta.copy();pt=tags==9
    overload=np.clip((e1[pt]-crit)/crit,0,5)
    new[pt]=np.clip(theta[pt]-.01*((theta[pt]-.5)/.5)**2*overload*.1,.5,1)
    assert np.all(new<=theta) and np.array_equal(new[~pt],theta[~pt])
    return new

class Model:
    def __init__(self,d,crit,pressure_mode):
        self.d=d;self.crit=crit;self.pressure_mode=pressure_mode
        assert d.MPI.size(d.MPI.comm_world)==1,'Serial execution only'
        report=json.loads((MESH/'conversion_report.json').read_text())
        for name in ('whole_eye_domain.xdmf','whole_eye_domain.h5','whole_eye_facet_region.xdmf',
                     'whole_eye_facet_region.h5','mesh_metadata.npz'):
            assert hashlib.sha256((MESH/name).read_bytes()).hexdigest()==report['files_sha256'][name]
        runtime_report=json.loads((RUNTIME_MESH/'verification.json').read_text())
        assert runtime_report['status']=='PASSED_MESH_AND_Q1_ROUNDTRIP'
        for name,digest in runtime_report['files_sha256'].items():
            assert hashlib.sha256((RUNTIME_MESH/name).read_bytes()).hexdigest()==digest
        mesh=d.Mesh()
        with d.XDMFFile(str(RUNTIME_MESH/'whole_eye_domain.xdmf')) as f:f.read(mesh)
        assert mesh.ufl_cell().cellname()=='hexahedron'
        assert mesh.num_cells()==36864 and mesh.num_vertices()==42494
        def read_tags(filename,dim):
            mvc=d.MeshValueCollection('size_t',mesh,dim)
            with d.XDMFFile(str(RUNTIME_MESH/filename)) as f:f.read(mvc,'name_to_read')
            return d.MeshFunction('size_t',mesh,mvc)
        self.tags=read_tags('whole_eye_domain.xdmf',3)
        raw=read_tags('whole_eye_facet_region.xdmf',2)
        self.facets=d.MeshFunction('size_t',mesh,2,0)
        valid=np.isin(raw.array(),[1,2,3,4,5]);self.facets.array()[valid]=raw.array()[valid]
        self.cell_tags=self.tags.array().copy();assert set(self.cell_tags)==set(MATERIALS)
        self.mesh=mesh;data=np.load(MESH/'mesh_metadata.npz')
        lookup={tuple(p):i for i,p in enumerate(data['points'])}
        vertex_map=np.array([lookup[tuple(p)] for p in mesh.coordinates()])
        key=lambda nodes:tuple(sorted(map(int,nodes)))
        cell_lookup={key(c):i for i,c in enumerate(data['hex_connectivity'])}
        cell_map=np.array([cell_lookup[key(vertex_map[c.entities(0)])] for c in d.cells(mesh)])
        assert len(set(cell_map))==mesh.num_cells()
        assert np.array_equal(self.cell_tags,data['material_id'][cell_map])
        inverse_map=np.empty(len(vertex_map),dtype=int);inverse_map[vertex_map]=np.arange(len(vertex_map))
        self.source_order_cells=inverse_map[data['hex_connectivity'][cell_map]]
        self.reference_nodes=mesh.coordinates()[self.source_order_cells]
        self.volume_cells=volumes(self.reference_nodes)
        mesh.init(2,3);mesh.init(2,0)
        face_lookup={key(c):int(t) for c,t in zip(data['boundary_quads'],data['boundary_id'])}
        seen=set();surface_vertices={t:set() for t in range(1,6)}
        for face in d.facets(mesh):
            k=key(vertex_map[face.entities(0)]);tag=int(self.facets[face])
            assert tag==face_lookup.get(k,0)
            if tag:
                assert face.exterior();seen.add(k);surface_vertices[tag].update(map(int,face.entities(0)))
        assert seen==set(face_lookup)
        self.V=d.VectorFunctionSpace(mesh,'Q',1)
        self.Q=d.FunctionSpace(mesh,'DQ',0)
        self.T=d.TensorFunctionSpace(mesh,'DQ',0,shape=(3,3))
        self.qmap=np.array([self.Q.dofmap().cell_dofs(i)[0] for i in range(mesh.num_cells())])
        self.tmap=np.array([[self.T.sub(j).dofmap().cell_dofs(i)[0] for j in range(9)] for i in range(mesh.num_cells())])
        assert len(np.unique(self.qmap))==self.Q.dim()
        assert len(np.unique(self.tmap))==self.T.dim()
        self.v2d=d.vertex_to_dof_map(self.V).reshape(-1,3)
        assert np.allclose(self.V.tabulate_dof_coordinates()[self.v2d],mesh.coordinates()[:,None,:])
        self.u=d.Function(self.V,name='Displacement')
        self.theta=d.Function(self.Q,name='Theta_Jg');self.set_cell(self.theta,np.ones(mesh.num_cells()))
        mu=d.Function(self.Q);lam=d.Function(self.Q)
        e=np.array([MATERIALS[int(t)][1] for t in self.cell_tags])
        self.set_cell(mu,e/(2*(1+NU)));self.set_cell(lam,e*NU/((1+NU)*(1-2*NU)))
        self.dx=d.Measure('dx',domain=mesh,subdomain_data=self.tags,metadata={'quadrature_degree':3})
        self.ds=d.Measure('ds',domain=mesh,subdomain_data=self.facets,metadata={'quadrature_degree':3})
        self.F=d.variable(d.Identity(3)+d.grad(self.u));self.J=d.det(self.F)
        Fg=self.theta**(1/3)*d.Identity(3);Fe=self.F*d.inv(Fg);Je=d.det(Fe)
        self.W=mu/2*(d.tr(Fe.T*Fe)-3)-mu*d.ln(Je)+lam/2*d.ln(Je)**2
        E=.5*(self.F.T*self.F-d.Identity(3));P=d.diff(self.W,self.F)
        sigma=P*self.F.T/self.J;dev=sigma-d.tr(sigma)/3*d.Identity(3)
        self.vm=d.sqrt(1.5*d.inner(dev,dev))
        self.factor=d.Constant(0);N=d.FacetNormal(mesh);v=d.TestFunction(self.V)
        # Default follows the previous kinematic branch: traction normal to REFERENCE surface.
        # Optional follower mode applies pressure normal to the CURRENT deformed surface.
        area_normal=N if pressure_mode=='reference' else self.J*d.inv(self.F).T*N
        self.R=d.derivative(self.W*self.dx,self.u,v)
        self.R+=self.factor*(.0047*d.dot(area_normal,v)*self.ds(3)+.00172*d.dot(area_normal,v)*self.ds(4))
        self.A=d.derivative(self.R,self.u,d.TrialFunction(self.V))
        self.bcs=[d.DirichletBC(self.V,d.Constant((0.,0.,0.)),self.facets,t) for t in (1,5)]
        self.bcs += [d.DirichletBC(self.V.sub(2),d.Constant(0.),self.facets,2)]
        actual=set().union(*(set(b.get_boundary_values()) for b in self.bcs))
        expected=set(self.v2d[sorted(surface_vertices[1]|surface_vertices[5])].ravel())
        expected.update(self.v2d[sorted(surface_vertices[2]),2])
        assert actual==expected;self.fixed=np.array(sorted(actual),dtype=int)
        self.Efield=d.Function(self.T)
        trial,test=d.TrialFunction(self.T),d.TestFunction(self.T)
        self.Esolver=d.LocalSolver(d.inner(trial,test)*self.dx,d.inner(E,test)*self.dx);self.Esolver.factorize()
        self.fields={};self.scalar_solvers={}
        for name,expr in [('J_total',self.J),('von_Mises',self.vm)]:
            f=d.Function(self.Q,name=name);a,b=d.TrialFunction(self.Q),d.TestFunction(self.Q)
            solver=d.LocalSolver(a*b*self.dx,expr*b*self.dx);solver.factorize()
            self.fields[name]=f;self.scalar_solvers[name]=solver
        for name in ('Material_ID','E1_total'):self.fields[name]=d.Function(self.Q,name=name)
        self.set_cell(self.fields['Material_ID'],self.cell_tags)
        measured={str(t):float(d.assemble(d.Constant(1)*self.dx(t))) for t in MATERIALS}
        for t in MATERIALS:
            assert np.isclose(measured[str(t)],report['material_volumes'][str(t)],rtol=1e-8,atol=1e-10)
        # Affine Q1 check verifies vector/component ordering and complete 3x3 projection.
        Ftest=np.array([[1.01,.02,.01],[.005,.99,.015],[0,.01,1.02]])
        uvalues=np.zeros(self.V.dim());uvalues[self.v2d]=mesh.coordinates()@(Ftest-np.eye(3)).T
        self.put(self.u,uvalues);tensor=self.tensor()
        tensor_error=float(np.max(abs(tensor-.5*(Ftest.T@Ftest-np.eye(3)))))
        assert tensor_error<1e-9
        self.put(self.u,np.zeros(self.V.dim()))
        assert d.assemble(self.R).norm('l2')<1e-8
        self.audit=dict(cells=mesh.num_cells(),vertices=mesh.num_vertices(),Q1_dofs=self.V.dim(),
            fixed_dofs=len(self.fixed),reference_volumes=measured,affine_tensor_error=tensor_error,
            material_counts={str(t):int(sum(self.cell_tags==t)) for t in MATERIALS},
            boundary_counts={str(t):int(sum(self.facets.array()==t)) for t in range(1,6)},
            pressure_mode=pressure_mode,quadrature_degree=3,cell_and_facet_mapping='PASS')

    def put(self,f,a):f.vector().set_local(a);f.vector().apply('insert')
    def set_cell(self,f,a):
        v=np.zeros(self.Q.dim());v[self.qmap]=a;self.put(f,v)
    def theta_cells(self):return self.theta.vector().get_local()[self.qmap]
    def tensor(self):
        self.Esolver.solve_local_rhs(self.Efield)
        a=self.Efield.vector().get_local()[self.tmap].reshape(-1,3,3)
        assert np.isfinite(a).all() and np.max(abs(a-a.transpose(0,2,1)))<1e-9
        return .5*(a+a.transpose(0,2,1))
    def e1(self):return np.linalg.eigvalsh(self.tensor())[:,2]
    def determinant_samples(self):
        nodes=self.mesh.coordinates()+self.u.vector().get_local()[self.v2d]
        return sampled_J(self.reference_nodes,nodes[self.source_order_cells])
    def residual(self):
        r=self.d.assemble(self.R)
        for b in self.bcs:b.apply(r)
        return r
    def solve(self):
        from bfgs_solver import solve
        return solve(self)

    def stats(self):
        d=self.d;theta=self.theta_cells()
        assert np.all(theta[self.cell_tags!=9]==1) and np.all((theta>=.5)&(theta<=1))
        assert np.max(abs(self.u.vector().get_local()[self.fixed]))<1e-12
        det=self.determinant_samples();assert np.isfinite(det).all() and det.min()>0
        data=dict(theta_min=float(theta.min()),sampled_J_min=float(det.min()),sampled_J_max=float(det.max()))
        data['sampled_Je_min']=float(np.min(det/theta[:,None]))
        for t,name in ((9,'PT'),(8,'LC')):
            data[name+'_volume']=float(d.assemble(self.J*self.dx(t)))
        pt=self.cell_tags==9
        data['PT_mean_theta']=float(self.volume_cells[pt]@theta[pt]/self.volume_cells[pt].sum())
        return data

    def write(self,file,index):
        for name,solver in self.scalar_solvers.items():solver.solve_local_rhs(self.fields[name])
        self.set_cell(self.fields['E1_total'],self.e1())
        for f in (self.u,self.theta,*self.fields.values()):file.write(f,float(index))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    modes=p.add_mutually_exclusive_group(required=True)
    modes.add_argument('--check-only',action='store_true');modes.add_argument('--smoke',action='store_true')
    modes.add_argument('--full-cycle',action='store_true')
    p.add_argument('--e1-threshold',type=float,default=.005)
    p.add_argument('--pressure-mode',choices=('reference',),default='reference')
    args=p.parse_args();assert np.isfinite(args.e1_threshold) and args.e1_threshold>0
    folder=ROOT/'results'/datetime.now().strftime('%Y%m%d_%H%M%S_%f');folder.mkdir(parents=True)
    print('Output:',folder,flush=True);(folder/'source').mkdir()
    for f in list(ROOT.glob('*.py'))+list(ROOT.glob('*.md')):shutil.copy2(f,folder/'source'/f.name)
    summary=dict(status='RUNNING_UNVERIFIED_EXPERIMENTAL',settings=vars(args),
                 materials={str(t):dict(name=n,E_MPa=e,nu=NU,placeholder=t in (4,5,6,7)) for t,(n,e) in MATERIALS.items()},
                 IOP_MPa=.0047,CSFP_MPa=.00172,states=[])
    summary['evolution']=dict(E1_crit=args.e1_threshold,theta_min=.5,growth_rate=.01,
                              gamma_theta=2.,stimulus_cap=5.,dt=.1,PT_tag=9)
    summary['source_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (folder/'source').iterdir()}
    summary['mesh_manifest']=json.loads((MESH/'conversion_report.json').read_text())
    summary['runtime_mesh_manifest']=json.loads((RUNTIME_MESH/'verification.json').read_text())
    save(folder/'summary.json',summary);files={};started=time.perf_counter()
    try:
        import dolfin as d
        assert d.has_lu_solver_method('mumps')
        m=Model(d,args.e1_threshold,args.pressure_mode);summary['setup']=m.audit
        with (folder/'measurements.csv').open('w',newline='') as logfile:
            writer=None
            def record(index,phase,step,solver):
                nonlocal writer
                row=dict(index=index,phase=phase,step=step,load_factor=float(m.factor),**m.stats())
                if writer is None:writer=csv.DictWriter(logfile,fieldnames=list(row));writer.writeheader()
                writer.writerow(row);logfile.flush()
                summary['states'].append(dict(**row,solver=solver));save(folder/'summary.json',summary)
            record(0,'reference',0,None)
            if not args.check_only:
                schedule=(10,3,10) if args.smoke else (100,500,100)
                index=0
                for phase,count in zip(('loading','remodeling','unloading'),schedule):
                    file=d.XDMFFile(str(folder/(phase+'.xdmf')));files[phase]=file
                    file.parameters['flush_output']=True;file.parameters['functions_share_mesh']=True
                    file.parameters['rewrite_function_mesh']=False
                    m.write(file,index);frozen=m.theta_cells().copy()
                    for step in range(1,count+1):
                        if phase=='loading':m.factor.assign(step/count)
                        elif phase=='unloading':m.factor.assign(1-step/count)
                        else:m.set_cell(m.theta,evolution(m.theta_cells(),m.e1(),m.cell_tags,args.e1_threshold))
                        solver=m.solve();index+=1
                        if phase!='remodeling':assert np.array_equal(m.theta_cells(),frozen)
                        record(index,phase,step,solver)
                        if step==1 or step%10==0 or step==count:m.write(file,index)
                        print(phase,step,count,solver,flush=True)
            summary['status']='SETUP_CHECK_PASSED' if args.check_only else 'COMPLETED_EXPLORATORY_NOT_VALIDATED'
    except Exception:
        summary['status']='FAILED';summary['traceback']=traceback.format_exc();raise
    finally:
        for file in files.values():file.close()
        summary['elapsed_seconds']=time.perf_counter()-started;save(folder/'summary.json',summary)
        print(summary['status'],folder,flush=True)

if __name__=='__main__':main()

