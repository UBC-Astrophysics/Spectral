from hydrogen_s import hydrogen_s
import numpy as np

betalist=np.linspace(0,1,101)
mlist=np.linspace(-5,5,11)
bb,mm = np.meshgrid(betalist,mlist)
bbf=bb.flatten()
mmf=mm.flatten()
Lamres=[]
for b,m in zip(bbf,mmf):
    V,Lam,w,rs,igood,zoom = hydrogen_s(b,11,31,m)
    Lamres.append(Lam)    
Lamres=np.reshape(Lamres,np.shape(bb))

linelist=[]
for i,b in enumerate(betalist):
    for j,m in enumerate(mlist)
        Lam1=Lamres[j,i]
        Lam1=Lam1[Lam1<0]
        Lam1=Lam1[:,np.newaxis]
        ll=[]
        for j2 in [j,j+1]:
            if (j2<len(mlist)):
                Lam2=Lamres[j2,i]
                Lam2=Lam2[Lam2<0]
                ll.append(np.flatten(Lam1-Lam2))
        linelist.append(ll)
linelist=np.reshape(linelist,np.shape(bb))

with open('linelist.dat','w') as f:
  for i,b in enumerate(betalist):
    for j,m in enumerate(mlist)
      f.write(b,m,linelist[i,j],'\n')
      
