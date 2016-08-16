% THREEDATOM - calculate an atom in three-D 
% the wavefunction is phi=f(r,mu=cos(theta))/r
% r coordinate, ranging from -1 to 1 
% Examples for Helium I
% nu=[1,1]; spin=[+1/2,-1/2]; mphi=[0, 0]; Z=2; %% 1s^2 1S  J=0 : 1.43189
% nu=[1,2]; spin=[-1/2,-1/2]; mphi=[0, 0]; Z=2; %% 1s2s 3S  J=1 : 1.08996
% nu=[1,2]; spin=[+1/2,-1/2]; mphi=[0, 0]; Z=2; %% 1s2s 1S  J=0 : 1.07930
% nu=[1,1]; spin=[-1/2,-1/2]; mphi=[0,-1]; Z=2; %% 1s2p 3P0 J=2 : 1.06856 
% nu=[1,3]; spin=[-1/2,-1/2]; mphi=[0, 0]; Z=2; %% 1s2p 3P0 J=1 : 1.06857 
% nu=[1,1]; spin=[+1/2,+1/2]; mphi=[0,-1]; Z=2; %% 1s2p 3P0 J=0 : 1.06856
% nu=[1,1]; spin=[-1/2,+1/2]; mphi=[0,-1]; Z=2; %% 1s2p 1P0 J=1 : 1.06593
% nu=[1,3]; spin=[-1/2,+1/2]; mphi=[0, 0]; Z=2; %% 1s2p 1P0 J=1 : 1.06594
% Examples for Lithium I
  nu=[1,1,2]; spin=[-1/2,+1/2,-1/2]; mphi=[0,0, 0]; Z=3; %% 1s^2 2s 2S  J=1/2 : 1.65595
% nu=[1,1,3]; spin=[-1/2,+1/2,-1/2]; mphi=[0,0, 0]; Z=3; %% 1s^2 2p 2P0 J=1/2 : 1.64097
% nu=[1,1,1]; spin=[-1/2,+1/2,-1/2]; mphi=[0,0,-1]; Z=3; %% 1s^2 2p 2P0 J=3/2 : 1.64097
% nu=[1,1,4]; spin=[-1/2,+1/2,-1/2]; mphi=[0,0, 0]; Z=3; %% 1s^2 3s 2S  J=1/2 : 1.62883
M=7; N=31;  MAXITER=100;  beta = 0;
% w coordinate, ranging from -1 to 1 
  [Dw,w] = cheb(M);  D2w = Dw^2;  Hw=-diag(1-w.*w)*D2w+diag(2.*w)*Dw;
% add zero weight to last term so we always use the same-sized vectors
  w_w=[inv(Dw(1:M,1:M))(1,:),0];
% r coordinate, ranging from 1 to -1, rp from 1 to 0, 
% rs from infinity to 0 (see below)
  [D,r] = cheb(N);  rp=0.5*(r+1);  D = 2*D;  D2 = D^2;
  hh=diag(1-rp.*rp);   D2=hh*(hh*D2+diag(-2*rp)*D);
  wr=inv(D(1:N,1:N))(1,:);
% drop infinite weight at infinity (i.e. 1) since w.f. goes to zero quickly
% now wr is the same dimension as D2 etc.
  wr=wr(2:N)./(1-rp(2:N).^2)';
% now multiply by hh so we didn't have zero for inversion above
  D=hh*D;
% Radial Laplacian minus r=infinity (for potential calcuations)
  D2inf = D2(2:N+1,2:N+1);
% Construct Laplacian minus r=0 and r=infinity 
% (set diricelet b.c. in radial direction for eigenvalue equations)
  D2 = D2(2:N,2:N); D =  D(2:N,2:N);  rp=rp(2:N);
  zoom=1/(1/116+sqrt(beta)/41)*0.5;  rs=zoom*atanh(rp);  wr=zoom*wr;
  R = diag(1./rs); R2 = diag(1./(rs.*rs));
% Hamilton for radial coordinate
  Hr=-1/(zoom*zoom)*D2-2*R;
  [rr,ww] = meshgrid(rs,w); rr = rr(:); ww = ww(:);  DIMEN=numel(rr);
  rperp2=rr.*rr.*(1-ww.*ww);  RR=diag(1./rr);  we=2*pi*kron(wr,w_w); 
% Full single electron Hamiltonian
  H=kron(Hr,eye(M+1))+kron(R2,Hw)+diag(beta*beta*rperp2);
  mmax=0;
  for e1cnt=1:numel(mphi)-1
    for e2cnt=e1cnt+1:numel(mphi)
      if (spin(e1cnt)==spin(e2cnt)) 
        dum=abs(mphi(e1cnt)-mphi(e2cnt));
        if (dum>mmax) mmax=dum; end
      end
    end
  end
% Construct Laplacians for electron potential calculation
  Linv=zeros(DIMEN,DIMEN,mmax+1);
  for i=0:mmax
% set the denominator equal to one, if it goes to zero (don't want an infinity)
% we will drop rperp2=0 later
    L=kron(1/(zoom*zoom)*D2,eye(M+1))-kron(R2,Hw)-diag(i*i./(rperp2+(rperp2==0)));
% set the BC for rperp2=0 but add in identity so can be inverted
% an infinity in the last step would have given an NaN here 
    if (i~=0) L=diag(rperp2==0)+diag(rperp2~=0)*L*diag(rperp2~=0); end
% add in the origin as a single point
    L=[ [L, kron(D2inf(1:N-1,N),ones(M+1,1))]; kron(D2inf(N,1:N-1),ones(1,M+1)), D2inf(N,N)];
% invert Laplacian and set boundary condition at origin
    Linvhold=(eye(DIMEN+1)-[zeros(DIMEN+1,DIMEN),ones(DIMEN+1,1)])*inv(L);
% divide by R before and after and remove origin!
% don't need to set B.C. again because the solutions for mphi~=0 are already zero in the 
% right places, and the e-val equation sets them to zero when needed
    Linv(:,:,i+1)=-4*pi*RR*Linvhold(1:DIMEN,1:DIMEN)*RR;
  end
  u=zeros(DIMEN,numel(nu));
  direct_old=zeros(DIMEN,DIMEN,numel(nu));
  exchange_old=zeros(DIMEN,DIMEN,numel(nu),mmax+1);
  direct_new=zeros(DIMEN,DIMEN,numel(nu));
  exchange_new=zeros(DIMEN,DIMEN,numel(nu),mmax+1);
  interact=zeros(DIMEN,DIMEN);
% begin iterations
  etotold=1e33;  etot=1;  j=1;
  while (abs((etot-etotold)/etot)>1e-5 & j<MAXITER) 
% calculate each electron wf
    etotold=etot;
    etot=0;
    for e1cnt=1:numel(nu)
      interact=zeros(DIMEN,DIMEN);
      for e2cnt=1:numel(nu)
        if (e2cnt~=e1cnt) 
          if (spin(e1cnt)==spin(e2cnt)) 
            interact=interact+direct_old(:,:,e2cnt) ...
                -exchange_old(:,:,e2cnt,abs(mphi(e1cnt)-mphi(e2cnt))+1);
          else
            interact=interact+direct_old(:,:,e2cnt);
          end
        end
      end
      if (mphi(e1cnt)==0) 
        [V,Lam] = eig(H+2/Z*interact);
      else
        Htmp=H+2/Z*interact+diag(mphi(e1cnt)*mphi(e1cnt)./rperp2);
% apply boundary condition along z-axis, remove columns and rows with rperp2==0
        [V,Lam] = eig(reshape(Htmp(find(rperp2*rperp2'~=0)),(N-1)*(M-1),(N-1)*(M-1)));
      end
      Lam = diag(Lam);  [Lam,ii] = sort(real(Lam)); 
      igood = find((V(1,:).*V(1,:))'<(M*N)^(-2)*1e-4);
      V = V(:,ii);   uu=V(:,igood(nu(e1cnt)));
% add in boundary condtion for angular direction (along z-axis), if needed
      if (mphi(e1cnt)~=0) 
        uu=[zeros(1,N-1);reshape(uu,M-1,N-1);zeros(1,N-1)](:); 
      end
      l(e1cnt)=Lam(igood(nu(e1cnt)))+2*beta*(mphi(e1cnt)+spin(e1cnt)-0.5);
      u2=uu.*uu; norm=we*u2; u2=u2/norm; uu=uu/sqrt(norm); u(:,e1cnt)=uu; % normalize
      eeen(e1cnt)=2/Z*we*(uu.*(interact*uu));  etot=etot+l(e1cnt)-eeen(e1cnt)/2;
      printf('J=%d e1=%d E-val: %18.12f EE: %18.12f Etot: %18.12f\n',j,e1cnt,l(e1cnt),eeen(e1cnt),etot);
      direct_new(:,:,e1cnt)=diag(Linv(:,:,1)*u2);
      for i=1:mmax+1 exchange_new(:,:,e1cnt,i)=diag(uu)*Linv(:,:,i)*diag(uu); end
    end
    j=j+1;  direct_old=direct_new;  exchange_old=exchange_new;
 end
 if (j==MAXITER) 
   printf('Did not converge\n');
 end

