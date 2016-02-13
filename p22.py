import numpy as np
from cheb import cheb
from scipy.linalg import eig
from scipy.special import airy
import matplotlib.pyplot as plt

# p22.py - 5th eigenvector of Airy equation u_xx = lambda*x*u

for N in range(12,60,12):
    D,x = cheb(N)
    D2 = np.dot(D,D)
    D2 = D2[1:-1,1:-1]
    Lam,V = eig(D2,np.diag(x[1:-1]))      # generalized ev problem
    Lam=np.real(Lam)
    V=np.real(V)
    V = V[:,Lam>0]
    Lam = Lam[Lam>0]
    iisort=np.argsort(Lam)
    ii=iisort[4]
    lam=Lam[ii]
    v = np.concatenate(([0],V[:,ii],[0]))
    v = v/v[N/2]*airy(0)[0]
    xx = np.linspace(-1,1,201)
    vv = np.polyval(np.polyfit(x,v,N),xx);
    plt.plot(xx,vv)
    plt.scatter(x,v)
    plt.title('N = %d     eig = %15.10f' % (N,lam))
    plt.xlim(-1,1)
    plt.show()
    plt.close()

