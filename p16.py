# p16.m - Poisson eq. on [-1,1]x[-1,1] with u=0 on boundary
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from cheb import cheb
from numpy import meshgrid, kron, eye, ravel,sin, dot, zeros, reshape, arange, linspace, vectorize
from numpy.linalg import solve
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

# Set up grids and tensor product Laplacian and solve for u:
N = 24
D,x = cheb(N)
y = x
xx,yy = meshgrid(x[1:-1],y[1:-1])
# stretch 2D grids to 1D vectors
xx = ravel(xx)
yy = ravel(yy)
f = 10*sin(8*xx*(yy-1));
D2 = dot(D,D)
D2 = D2[1:-1,1:-1]
I = eye(N-1)
L = kron(I,D2) + kron(D2,I)     # Laplacian
u=solve(L,f)

# Reshape long 1D results onto 2D grid:
uu = zeros((N+1,N+1))
ii=arange(1,N)
uu[meshgrid(ii,ii)] = reshape(u,(N-1,N-1))
xx,yy = meshgrid(x,y)
value = uu[N/4+1,N/4+1]

# Interpolate to finer grid and plot:
xxx=linspace(-1,1,51)
[xxx,yyy] = meshgrid(xxx,xxx)
ifunk=vectorize(interp2d(xx,yy,uu))
uuu=ifunk(xxx,yyy)  
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xxx,yyy,uuu, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()
