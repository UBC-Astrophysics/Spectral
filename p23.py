import numpy as np
from cheb import cheb
from scipy.linalg import eig
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import math

# p23.m - eigenvalues of perturbed Laplacian on [-1,1]x[-1,1]
#         (compare p16.m)

# Set up tensor product Laplacian and compute 4 eigenmodes:
N = 16
D,x = cheb(N)
y = x
xx,yy = np.meshgrid(x[1:-1],y[1:-1])
xx = np.ravel(xx)
yy = np.ravel(yy)
D2 = np.dot(D,D)
D2 = D2[1:-1,1:-1]
I = np.eye(N-1)
L = -np.kron(I,D2) - np.kron(D2,I)             #  Laplacian
L = L + np.diag(np.exp(20*(yy-xx-1)))          #  + perturbation
D,V = eig(L)
D=np.real(D)
V=np.real(V)
isort=np.argsort(D)
V=V[:,isort]
D=D[isort]

# interpolate to finer grid, and plot:
fine = np.linspace(-1,1,101)

# apply bc of zero along edges
uu = np.zeros((N+1,N+1))
for i in range(0,4):
  # set interior values
  uu[1:-1,1:-1] = np.reshape(V[:,i],(N-1,N-1))

  uu = uu/np.max(np.abs(uu))
  funk = interp2d(x,y,uu,kind='cubic');
  uuu=np.reshape(funk(fine,fine),(len(fine),len(fine)))
  
  plt.contour(fine,fine,uuu)
  plt.title('eig = %18.14f pi^2/4' % (D[i]/(math.pi**2/4)))
  plt.show()
  plt.close()
#     subplot('position',[ax(i) ay(i) .38 .38])
#     contour(fine,fine,uuu,-.9:.2:.9)
#     colormap(1e-6*[1 1 1]); axis square

#   end

