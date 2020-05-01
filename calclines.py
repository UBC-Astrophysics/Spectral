from hydrogen_s import hydrogen_s
import numpy as np
import matplotlib.pyplot as plt

betalist=np.linspace(0,1,101)
mlist=np.linspace(-6,6,13)
bb,mm = np.meshgrid(betalist,mlist)
bbf=bb.flatten()
mmf=mm.flatten()
Lamres=[]
for b,m in zip(bbf,mmf):
    V,Lam,w,rs,igood,zoom = hydrogen_s(b,11,31,m)
    Lam=np.real(Lam)
    Lam=-np.exp(np.unique(np.round(np.log(-Lam[Lam<-1.0/36.0]),decimals=3)))
    Lamres.append(Lam)    

Lamres=np.reshape(Lamres,np.shape(bb))
for i,b in enumerate(betalist):
    for j,m in enumerate(mlist):
        Lam1=np.real(Lamres[j,i])
        Lam1=Lam1[Lam1<-1.0/36.0]
        Lam1=Lam1[:,np.newaxis]
        for j2 in [j-1,j,j+1]:
            if (j2<len(mlist)) and (j2>=0):
                Lam2=np.real(Lamres[j2,i])
                Lam2=Lam2[Lam2<-1.0/36.0]
                res=(Lam1-Lam2).flatten()
                res=res[res>0.09]
                #res=res[(res>0.09) & (res<0.3)]
                res=np.exp(-np.unique(-np.round(np.log(res[res>0]),decimals=4)))
                print(('%g %d %d'+' %g'*len(res)) % ((b,m,mlist[j2])+tuple(res)))
                plt.plot(9.1164851/res,b+0*res,'b,')
plt.savefig('lines.pdf')
