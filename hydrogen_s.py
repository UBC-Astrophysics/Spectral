# HYDROGEN_S - calculate hydrogen in spherical coordinates
# beta is B/Bc, M is number of mu, N is number of radial
# the wavefunction is phi=f(r, mu=cos(theta))/r exp(i m phi)
# eigenfunctions in V, e-val in Lam, coords in (w,rs)

from cheb import cheb
import numpy as np
from scipy.linalg import eig

def hydrogen_s(beta,M,N,mphi):
# w coordinate, ranging from -1 to 1
    Dw,w = cheb(M)
    D2w = np.dot(Dw,Dw)
    if (mphi!=0):
        w=w[1:-1]
        Dw=Dw[1:-1,1:-1]
        D2w=D2w[1:-1,1:-1]
    Hw=-np.dot(np.diag(1-w*w),D2w)+np.dot(np.diag(2.*w),Dw);
# r coordinate, ranging from 1 to -1, rp from 1 to 0
    D,r = cheb(N)
    rp = 0.5*(r+1)
    D = 2*D
    D2 = np.dot(D,D)
    hh = np.diag(1-rp*rp)
    D2 = np.dot(hh,( np.dot(hh,D2)+ np.dot(np.diag(-2*r),D)))
    D = np.dot(hh,D)
    D2 = D2[1:-1,1:-1]
    D = D[1:-1,1:-1]
    rp = rp[1:-1]
# zoom factor: set by coulomb and larmor radius; , rs from inf to 0
    zoom=1/(1.0/110.0+beta**0.5/41)
    rs=zoom*np.arctanh(rp)
    R = np.diag(1/rs)
    R2 = np.diag(1/(rs*rs))
    Hr=-1/(zoom*zoom)*D2-2*R
    rr,ww = np.meshgrid(rs,w)
    rr = np.ravel(rr)
    ww = np.ravel(ww)
    rperp2=rr*rr*(1-ww*ww)
    if (mphi==0):
        H = np.kron(Hr,np.eye(M+1))+np.kron(R2,Hw)+np.diag(beta*beta*rperp2)
    else:
        H = np.kron(Hr,np.eye(M-1))+np.kron(R2,Hw)+np.diag(beta*beta*rperp2)+np.diag(mphi*mphi/rperp2)
    Lam,V = eig(H)
    ii=np.argsort(Lam)
    Lam = Lam[ii]
    Lam = Lam+2*beta*(mphi-1)
    V = V[ii]
#     check outer B.C. and for bound states
    igood=0
#  igood = find((V(1,:).*V(1,:))'<(M*N)^(-2)*1e-4 & Lam<0);
# 
    return V,Lam,w,rs,igood,zoom
