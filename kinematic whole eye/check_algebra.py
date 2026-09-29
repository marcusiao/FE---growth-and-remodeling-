"""No tissue simulation: BFGS algebra and retained DOLFIN LU checks."""
import json
import numpy as np
import dolfin as d
from petsc4py import PETSc
PETSc.Log.begin()
from bfgs_solver import inverse_action

rng=np.random.default_rng(184)
B=rng.normal(size=(12,12)); B=B.T@B+np.eye(12)
H=np.linalg.inv(B); pairs=[]; errors=[]
for i in range(10):
    s=rng.normal(size=12);y=B@s+.1*s;rho=1/float(s@y)
    pairs.append((s,y,rho));V=np.eye(12)-rho*np.outer(s,y)
    H=V@H@V.T+rho*np.outer(s,s)
    v=rng.normal(size=12)
    errors.append(float(np.linalg.norm(inverse_action(v,lambda q:np.linalg.solve(B,q),pairs)-H@v)))
assert max(errors)<1e-11
mesh=d.UnitIntervalMesh(8);Q=d.FunctionSpace(mesh,'CG',1)
u=d.TrialFunction(Q);v=d.TestFunction(Q)
A=d.assemble((u*v+d.inner(d.grad(u),d.grad(v)))*d.dx)
lu=d.LUSolver(A,'mumps')
res=[]
for scale in (1.,2.,-.3):
    rhs=d.assemble(d.Constant(scale)*v*d.dx);x=rhs.copy();lu.solve(x,rhs)
    e=rhs.copy();A.mult(x,e);e.axpy(-1,rhs);res.append(e.norm('l2'))
assert max(res)<1e-11
factor_count=PETSc.Log.Event('MatLUFactorNum').getPerfInfo()['count']
assert factor_count==1, ('Expected one numeric factorization for three solves',factor_count)
report=dict(status='PASSED_ALGEBRA_ONLY_NO_MECHANICAL_CYCLE',max_BFGS_error=max(errors),
            reused_LU_residuals=res,numeric_factorizations=factor_count)
from pathlib import Path
Path(__file__).with_name('algebra_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report))

