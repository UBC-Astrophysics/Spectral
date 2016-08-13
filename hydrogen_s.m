% HYDROGEN_S - calculate hydrogen in spherical coordinates
% beta is B/Bc, M is number of mu, N is number of radial
% the wavefunction is phi=f(r, mu=cos(theta))/r exp(i m phi)
% eigenfunctions in V, e-val in Lam, coords in (w,rs)
function [V,Lam,w,rs,igood,zoom] = hydrogen_s(beta,M,N,mphi)
% w coordinate, ranging from -1 to 1
  [Dw,w] = cheb(M); D2w = Dw^2;
  if (mphi~=0) w=w(2:M); Dw=Dw(2:M,2:M); D2w=D2w(2:M,2:M); end
  Hw=-diag(1-w.*w)*D2w+diag(2.*w)*Dw;
% r coordinate, ranging from 1 to -1, rp from 1 to 0
  [D,r] = cheb(N); rp = 0.5*(r+1); D = 2*D; D2 = D^2;
  hh = diag(1-rp.*rp); D2 = hh*(hh*D2+diag(-2*r)*D); D = hh*D;
  D2 = D2(2:N,2:N); D = D(2:N,2:N); rp = rp(2:N);
% zoom factor: set by coulomb and larmor radius; , rs from inf to 0
  zoom=1/(1/110+sqrt(beta)/41); rs=zoom*atanh(rp); R = diag(1./rs); R2 = diag(1./(rs.*rs));
  Hr=-1/(zoom*zoom)*D2-2*R; [rr,ww] = meshgrid(rs,w); rr = rr(:); ww = ww(:); rperp2=rr.*rr.*(1-ww.*ww);
  if (mphi==0)
    H = kron(Hr,eye(M+1))+kron(R2,Hw)+diag(beta*beta*rperp2);
  else
    H = kron(Hr,eye(M-1))+kron(R2,Hw)+diag(beta*beta*rperp2)+diag(mphi*mphi./rperp2);
  end
  [V,Lam] = eig(H); Lam = diag(Lam); [Lam,ii] = sort(Lam); Lam = Lam+2*beta*(mphi-1); V = V(:,ii);
% check outer B.C. and for bound states
  igood = find((V(1,:).*V(1,:))'<(M*N)^(-2)*1e-4 & Lam<0);
