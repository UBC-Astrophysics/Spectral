from hydrogen_s import hydrogen_s
import numpy as np
import matplotlib.pyplot as plt

betalist=np.linspace(0,0.1,201)
mmax=5
mlist=np.linspace(-mmax,mmax,2*mmax+1)
mm,bb = np.meshgrid(mlist,betalist)
bbf=bb.flatten()
mmf=mm.flatten()
Lamres=[]
for b,m in zip(bbf,mmf):
    V,Lam,w,rs,igood,zoom = hydrogen_s(b,11,31,m)
    Lam=np.sort(np.real(Lam))
    Lam=Lam[0:(mmax-(abs(m))+2)*(mmax-abs(m)+1)/2]
    if (len(Lam)>0):
        print(('# %g %d %d'+' %g'*len(Lam)) % ((b,m,99)+tuple(Lam)))
    Lamres.append(Lam)    

Lamres=np.reshape(Lamres,np.shape(bb))
for i,b in enumerate(betalist):
    resbeta=[]
    for j,m in enumerate(mlist):
        if (m==0) or (m==-1) or (m==1):
            Lam1=Lamres[i,j]
            Lam1=Lam1[:,np.newaxis]
            for j2 in [j-1,j,j+1]:
                if (j2<len(mlist)) and (j2>=0):
                    Lam2=Lamres[i,j2]
                    res=(Lam1-Lam2).flatten()
                    resbeta=np.concatenate((resbeta,res))
                    if (len(res)>0):
                        print(('# %g %d %d'+' %g'*len(res)) % ((b,m,mlist[j2])+tuple(res)))
    print(('%g %d %d'+' %g'*len(resbeta)) % ((b,-99,-99)+tuple(resbeta)))
