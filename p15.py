import matplotlib.pyplot as plt
from cheb import cheb
import numpy as np
from scipy.linalg import eig
import math as mt

# p15.m - solve eigenvalue BVP u_xx = lambda*u, u(-1)=u(1)=0
N = 36
D,x = cheb(N)
D2 = np.dot(D,D)
D2 = D2[1:-1,1:-1]
lam,V = eig(D2)
lam=np.real(lam)
V=np.real(V)
ii=np.argsort(-lam)
V = V[:,ii]
lam = lam[ii]
eiglist=range(4,34,5)
for p,j in enumerate(eiglist):
    plt.subplot(len(eiglist),1,p+1)
    u = np.concatenate(([0],V[:,j],[0]))
    plt.scatter(x,u)
    xx = np.linspace(-1,1,201)
    uu = np.polyval(np.polyfit(x,u,N),xx)
    plt.plot(xx,uu)
    plt.title('eig %d =%20.13f*pi^2/4' % (j+1,lam[j]*4/mt.pi**2))

plt.show()

