import glob,re,sys,numpy as np
def arr(t,n,dt=float):
    m=re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>'%n,t,re.S); return np.array(m.group(1).split(),dtype=dt)
num=den=0.0; zs=[]
for f in sorted(glob.glob(sys.argv[1])):
    t=open(f,errors='replace').read(); P=arr(t,'Points').reshape(-1,3); c=arr(t,'connectivity',int).reshape(-1,8); u=arr(t,'Velocity').reshape(-1,3)
    X=P[c]; vol=np.zeros(len(c))
    for i,j,k,l in [(0,1,3,4),(1,2,3,6),(1,3,4,6),(1,4,5,6),(3,4,6,7)]:
        vol+=np.abs(np.einsum('ij,ij->i',np.cross(X[:,j]-X[:,i],X[:,k]-X[:,i]),X[:,l]-X[:,i]))/6
    num+=np.sum(u[c][:,:,2].mean(axis=1)*vol); den+=vol.sum()
print(sys.argv[1], "cells-volume", den, "unweighted mean u_z", num/den)
