import numpy as np
from cheb import cheb
from numpy.linalg import norm
import matplotlib.pyplot as plt

# p12.py - accuracy of Chebyshev spectral differentiation
#         (compare p7.py)

# Compute derivatives for various values of N:
Nmax = 50
E=np.zeros((5,Nmax))
for N in range(1,Nmax+1):
    D,x = cheb(N)
    v = np.abs(x)**3
    vprime = 3*x*np.abs(x)                   # 3rd deriv in BV
    E[1,N-1] = norm(np.dot(D,v)-vprime,np.inf);
    v = np.exp(-1/(x*x))
    vprime = 2*v/(x*x*x)                     # C-infinity
    E[2,N-1] = norm(np.dot(D,v)-vprime,np.inf);
    v = 1/(1+x*x)
    vprime = -2*x*v*v                        # analytic in [-1,1]
    E[3,N-1] = norm(np.dot(D,v)-vprime,np.inf);
    v = x**10
    vprime = 10*x**9                         # polynomial
    E[4,N-1] = norm(np.dot(D,v)-vprime,np.inf);

# Plot results:
titles = ['|x^3|','exp(-x^{-2})','1/(1+x^2)','x^{10}']
for iplot in range(1,5):
    plt.semilogy(np.arange(1,Nmax+1),E[iplot])
    plt.title(titles[iplot-1])
    plt.show()
    plt.close()
