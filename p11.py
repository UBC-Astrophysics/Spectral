from cheb import cheb
import numpy as np
import matplotlib.pyplot as plt

# p11.py - Chebyshev differentation of a smooth function

xx = np.linspace(-1,1,201)
uu = np.exp(xx)*np.sin(5*xx)
for N in (10,20):
    D,x = cheb(N)
    u = np.exp(x)*np.sin(5*x)
    plt.plot(xx,uu)
    plt.scatter(x,u)
    plt.xlim(-1,1)
    plt.show()
    plt.close()
    duu=np.exp(xx)*(np.sin(5*xx)+5*np.cos(5*xx))
    plt.plot(xx,duu)
    du=np.dot(D,u)
    plt.scatter(x,du)
    plt.xlim(-1,1)
    plt.show()
    plt.close()
    error = du-np.exp(x)*(np.sin(5*x)+5*np.cos(5*x))
    plt.scatter(x,error)
    plt.show()
    plt.xlim(-1,1)
    plt.close()
