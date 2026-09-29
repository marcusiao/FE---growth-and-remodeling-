"""Inverse BFGS with a retained exact-tangent LU base; reference loads only.

FEBio-style reformation/reuse strategy, not a port of FEBio's solver.
Retains the baseline force-residual line search and convergence criteria.
"""
import time
import numpy as np


def inverse_action(vector, base_solve, pairs):
    q=vector.copy(); coefficients=[]
    for s,y,rho in reversed(pairs):
        a=rho*np.dot(s,q); coefficients.append(a); q-=a*y
    z=base_solve(q)
    for (s,y,rho),a in zip(pairs,reversed(coefficients)):
        z+=s*(a-rho*np.dot(y,z))
    return z


def solve(model):
    if model.pressure_mode!='reference':
        raise ValueError('BFGS candidate supports symmetric reference-load problem only')
    d=model.d; start=time.perf_counter()
    timing=dict(determinant=0.,residual=0.,assembly=0.,base_linear_solve=0.,bfgs_action=0.)
    counts=dict(reformations=0,base_solves=0,updates=0,line_search_trials=0)
    history=[]; events=[]; pairs=[]; lu=None; A=None
    def timed(key,fn):
        t=time.perf_counter(); value=fn();timing[key]+=time.perf_counter()-t;return value
    def valid():
        j=timed('determinant',model.determinant_samples)
        return np.isfinite(j).all() and j.min()>1e-8
    def residual():return timed('residual',model.residual)
    def reform(reason):
        nonlocal lu,A
        if counts['reformations']>=15:raise RuntimeError('BFGS reformation limit')
        A=timed('assembly',lambda:d.assemble(model.A))
        for bc in model.bcs:bc.apply(A)
        lu=d.LUSolver(A,'mumps')
        pairs.clear();counts['reformations']+=1;events.append(reason)
    def base_solve(values):
        rhs=model.u.vector().copy();rhs.set_local(values);rhs.apply('insert')
        out=rhs.copy();out.zero()
        timed('base_linear_solve',lambda:lu.solve(out,rhs));counts['base_solves']+=1
        error=rhs.copy();A.mult(out,error);error.axpy(-1,rhs)
        if error.norm('l2')/max(rhs.norm('l2'),1e-30)>1e-6:
            raise RuntimeError('Base LU residual check failed')
        result=out.get_local();result[model.fixed]=0.;return result
    for bc in model.bcs:bc.apply(model.u.vector())
    if not valid():raise RuntimeError('Invalid initial sampled J')
    r=residual();initial=r.norm('l2');tol=max(1e-9,1e-8*initial)
    for it in range(31):
        err=r.norm('l2');history.append(err)
        if np.isfinite(err) and err<=tol:
            return dict(iterations=it,residual=err,tolerance=tol,seconds=time.perf_counter()-start,
                        timing=timing,counts=counts,reformation_reasons=events,residual_history=history,
                        method='inverse_BFGS_exact_LU_base_max10')
        if not np.isfinite(err) or it==30:raise RuntimeError('BFGS convergence failure')
        if lu is None:reform('initial or curvature restart')
        elif len(pairs)>=10:reform('10 updates reached')
        old=model.u.vector().get_local();rv=r.get_local();rv[model.fixed]=0.
        accepted=False
        for attempt in range(2):
            delta=-timed('bfgs_action',lambda:inverse_action(rv,base_solve,pairs))
            delta[model.fixed]=0.
            for k in range(24):
                counts['line_search_trials']+=1
                model.put(model.u,old+(.5**k)*delta)
                if valid():
                    trial=residual();trialnorm=trial.norm('l2')
                    if np.isfinite(trialnorm) and trialnorm<=(1-1e-4*(.5**k))*err:
                        accepted=True;break
            if accepted:break
            model.put(model.u,old)
            if not pairs:break
            reform('updated-direction line search failed')
        if not accepted:
            model.put(model.u,old);raise RuntimeError('BFGS line search failed')
        s=model.u.vector().get_local()-old;y=trial.get_local()-rv
        s[model.fixed]=0.;y[model.fixed]=0.;sy=float(np.dot(s,y))
        if np.isfinite(sy) and sy>1e-12*np.linalg.norm(s)*np.linalg.norm(y):
            pairs.append((s,y,1/sy));counts['updates']+=1
        else:
            pairs.clear();lu=None;events.append('curvature safeguard requested restart')
        r=trial

