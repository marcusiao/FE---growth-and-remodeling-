"""NumPy trilinear hex8 geometry in meshio/VTK source node order."""
import numpy as np

SIGNS=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],dtype=float)

def derivatives(points):
    p=np.asarray(points)
    out=np.empty((len(p),8,3))
    for j in range(3):
        other=[k for k in range(3) if k!=j]
        out[:,:,j]=SIGNS[None,:,j]/8*np.prod(1+p[:,None,other]*SIGNS[None,:,other],axis=2)
    return out

def gauss(order=3):
    x,w=np.polynomial.legendre.leggauss(order)
    points=np.array([(a,b,c) for a in x for b in x for c in x])
    weights=np.array([a*b*c for a in w for b in w for c in w])
    return points,weights

def jacobians(nodes,points):
    return np.einsum('cai,qaj->cqij',nodes,derivatives(points))

def volumes(nodes):
    points,weights=gauss()
    det=np.linalg.det(jacobians(nodes,points))
    if not np.all(det>0): raise ValueError('Nonpositive reference hex Jacobian')
    return det@weights

def sampled_J(reference_nodes,current_nodes):
    points=np.vstack((gauss()[0],SIGNS,np.zeros((1,3))))
    a=np.linalg.det(jacobians(reference_nodes,points))
    b=np.linalg.det(jacobians(current_nodes,points))
    if not np.all(a>0): raise ValueError('Invalid reference hex geometry')
    return b/a
