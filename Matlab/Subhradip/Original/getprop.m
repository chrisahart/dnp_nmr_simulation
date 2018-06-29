function[prop]=getprop(evecham,evalgham,h_mw,spinops,t1e,t2e,t1n,t2n,we,wn,wr,nsteps)

S1z=spinops.S1z;
S2z=spinops.S2z;
Iz=spinops.Iz;
S1p=spinops.S1p;
S2p=spinops.S2p;
Ip=spinops.Ip;
S1m=spinops.S1m;
S2m=spinops.S2m;
Im=spinops.Im;
T=100;
trstep=(1/wr)/nsteps;
nsp=log(length(S1z))/log(2);
kb=1.3806e-23;
pl=planck;
t_corr=pl/(kb*T);
p_e=tanh(.5*we*t_corr);
p_n=tanh(.5*wn*t_corr);
gnp=0.5*(1-p_n)*(1/t1n);
gnm=0.5*(1+p_n)*(1/t1n);
gep=0.5*(1-p_e)*(1/t1e);
gem=0.5*(1+p_e)*(1/t1e);


prop=zeros(4^nsp,4^nsp,nsteps);
parfor ii=1:nsteps
    D=squeeze(evecham(:,:,ii));
    Dinv=D\eye(2^nsp); %%%% same as inv(D)
    hamil_mwt=Dinv*h_mw*D; %%% Transform the MW Hamiltonian to the same frame
    hamil_0=diag(evalgham(:,ii));
    hamilt=hamil_0+hamil_mwt; %%%% Total Hamiltonian
 
 %%%%%% Tilt all the operators to the same frame %%%%

    S1zt=Dinv*S1z*D;
    Izt=Dinv*Iz*D;
    S2zt=Dinv*S2z*D;
    S1pt=Dinv*S1p*D;
    S2pt=Dinv*S2p*D;
    S1mt=Dinv*S1m*D;
    S2mt=Dinv*S2m*D;
    Ipt=Dinv*Ip*D;
    Imt=Dinv*Im*D;
    
    
   %%%%%% Make left and right Liouvillian %%%%%%%%% Generate the
   %%%%%% Lindbladians
    LS1zt=kron(S1zt,S1zt.');
    RS1zt=kron(eye(2^nsp),eye(2^nsp));
    LS2zt=kron(S2zt,S2zt.');
    RS2zt=kron(eye(2^nsp),eye(2^nsp));
    LIzt=kron(Izt,Izt.');
    RIzt=kron(eye(2^nsp),eye(2^nsp));
    LIpt=kron(Ipt,Imt.')-.5*eye(4^nsp)+.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    LImt=kron(Imt,Ipt.')-.5*eye(4^nsp)-.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    LS1pt=kron(S1pt,S1mt.')-.5*eye(4^nsp)+.5*(kron(S1zt,eye(2^nsp))+kron(eye(2^nsp),S1zt.'));
    LS1mt=kron(S1mt,S1pt.')-.5*eye(4^nsp)-.5*(kron(S1zt,eye(2^nsp))+kron(eye(2^nsp),S1zt.'));
    LS2pt=kron(S2pt,S2mt.')-.5*eye(4^nsp)+.5*(kron(S2zt,eye(2^nsp))+kron(eye(2^nsp),S2zt.'));
    LS2mt=kron(S2mt,S2pt.')-.5*eye(4^nsp)-.5*(kron(S2zt,eye(2^nsp))+kron(eye(2^nsp),S2zt.'));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Rt2ematv3=(2/t2e)*(LS1zt-.25*RS1zt);
    Rt2e2matv3=(2/t2e)*(LS2zt-.25*RS2zt);   
    Rt2nmatv3=(2/t2n)*(LIzt-.25*RIzt);
    Rt1mat=gep*LS1pt+gem*LS1mt+gep*LS2pt+gem*LS2mt+...               %(-2*LSyt*RSyt+LSyt*LSyt+RSyt*RSyt)+(-2*LSzt*RSzt+LSzt*LSzt+RSzt*RSzt))+...
        gnp*LIpt+gnm*LImt;
    Rtot=1*Rt1mat+1*Rt2ematv3+1*Rt2e2matv3+1*Rt2nmatv3;

    Lhamilt=kron((hamilt),eye(2^nsp))-kron(eye(2^nsp),(hamilt).');
 
    LD=kron(D,D);%-kron(eye(4),D.');
    LDinv=kron(Dinv,Dinv);%-kron(eye(4),Dinv.');
    Ltot=Lhamilt+1*1i*((Rtot));
    prop(:,:,ii)=LD*expm(-1i*Ltot*trstep)*LDinv;
end
