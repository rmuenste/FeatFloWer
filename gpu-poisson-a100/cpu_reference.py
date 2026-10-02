#!/usr/bin/env python3
"""CPU reference for a synthetic instance: spectrum bounds (inverse iteration + Lanczos),
Jacobi-CG and 4x4 block-Jacobi-CG iteration counts to 1e-8 on a rough mean-free RHS, and
the solve of the stored b against x_true. Writes results/<instance>.cpu.json."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
f=sys.argv[1]; d=np.load(f); A=sp.csr_matrix((d["data"],d["indices"],d["indptr"])); nu=A.shape[0]; h=float(d["h"])
Dinv=sp.diags(1/A.diagonal())
x=np.random.default_rng(0).standard_normal(nu); x/=np.linalg.norm(x)
for it in range(10):
    y,_=spla.cg(A,x,rtol=1e-9,maxiter=50000,M=Dinv); x=y/np.linalg.norm(y)
lmin=float(x@(A@x)); lmax=float(spla.eigsh(A,k=1,which="LA",return_eigenvectors=False,tol=1e-4)[0])
one=np.ones(nu); rng=np.random.default_rng(1); b=rng.standard_normal(nu); b-=one*(one@b)/nu
out={"instance":f,"h":h,"nu":nu,"nnz":int(A.nnz),"nnz_per_row":A.nnz/nu,"lambda_min":lmin,"lambda_max":lmax,"condition":lmax/lmin,"h^-2":1/h/h}
def run(name,M):
    n=[0]; t=time.time(); xs,info=spla.cg(A,b,rtol=1e-8,maxiter=100000,M=M,callback=lambda xk:n.__setitem__(0,n[0]+1))
    out[name]={"iters":n[0],"info":int(info),"seconds":time.time()-t,"final_relres":float(np.linalg.norm(b-A@xs)/np.linalg.norm(b))}
    print(name, out[name], flush=True)
run("cg_jacobi",Dinv)
nc=nu//4; blocks=[sp.csr_matrix(np.linalg.inv(A[4*c:4*c+4,4*c:4*c+4].toarray())) for c in range(nc)]
run("cg_block_jacobi_4x4", sp.block_diag(blocks).tocsr())
xs,_=spla.cg(A,d["b"],rtol=1e-10,maxiter=100000,M=Dinv); out["smooth_rhs_rel_err_vs_x_true"]=float(np.linalg.norm(xs-d["x_true"])/np.linalg.norm(d["x_true"]))
print(json.dumps(out,indent=1)); json.dump(out,open("results/"+f.split("/")[-1].replace(".npz",".cpu.json"),"w"),indent=1)
