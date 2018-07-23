function[iznew,sznew,ps0,pi0]=fn_se_buildup(a,b,g,nr,final_time,mw_pwr,wme,t1e,t1n)
%tic
spins = [1/2 1/2];
nsp = length(spins);
Sz=sop(spins,'ze'); % You need easyspin sop file. 
Iz=sop(spins,'ez');
Ix=sop(spins,'ex');
Iy=sop(spins,'ey');
IzSz=sop(spins,'zz');
IpSz=sop(spins,'z+');
ImSz=sop(spins,'z-');
Sp=sop(spins,'+e');
Sm=sop(spins,'-e');
Ip=sop(spins,'e+');
Im=sop(spins,'e-');
Sx=sop(spins,'xe');
Sy=sop(spins,'ye');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Define the constants and magnitudes %%%%%%%
wn=-400.9e6;
we=28.025*9.4*1e9;
%we=263e9;
%we=wme-wn;
w1=mw_pwr*1e6;
delw=we-wme;
spinning_freq=nr;
wr = spinning_freq*1000;
tr = 1/wr;
nsteps = 10000;
trstep = tr/nsteps;
tarray=0:trstep:tr;
hzz_max=3e6;
%rho0 = Sz;
beta_en=b;
gamma_en=g;
h_mw=w1*Sx;
hbar=1.054e-34;
planckc=hbar*2*pi;
kb=1.3806e-23;
T=100;
% t1e=0.3e-3;
% t1n=10;
t2e=1e-6;
t2n=1e-3;
pl=planck;
t_corr=pl/(kb*T);
p_e=tanh(0.5*we*t_corr); 
p_n=tanh(0.5*wn*t_corr);
gnp=0.5*(1-p_n)*(1/(1*t1n)); %%%% This is thermalisation, expansion of exp-(t_corr*H)/Trace(same);
gnm=0.5*(1+p_n)*(1/(1*t1n));
gep=0.5*(1-p_e)*(1/(1*t1e));
gem=0.5*(1+p_e)*(1/(1*t1e));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% define g-anisotropy %%%%%%%%
thetam=acosd(1/sqrt(3));
alpha_g=a;
beta_g=b;
gamma_g=g;
gx=we*(2.00614/2)+18.76e6;
gy=we*(2.00194/2)+92.4e6;
gz=we*(2.00988/2)+18.2e6;
delta=gz-we;
eta=(gy-gx)/delta;

ca = cosd(alpha_g);
cb = cosd(beta_g);
cg = cosd(gamma_g);
sa = sind(alpha_g);
sb = sind(beta_g);
sg = sind(gamma_g);

r11 = ca*cb*cg-sa*sg;
r12 = sa*cb*cg+ca*sg; 
r13 = -sb*cg;
r21 = -ca*cb*sg-sa*cg; 
r22 = -sa*cb*sg+ca*cg; 
r23 = sb*sg;
r31 = ca*sb; 
r32 = sa*sb;
r33 = cb;
c0 = 1/3*(gx+gy+gz);
c1 = 2*sqrt(2)/3*(gx*r11*r31+gy*r12*r32+gz*r13*r33);
c2 = 2*sqrt(2)/3*(gx*r21*r31+gy*r22*r32+gz*r23*r33);
c3 = 1/3*(gx*(r11^2-r21^2)+gy*(r12^2-r22^2)+gz*(r13^2-r23^2));
c4 = 2/3*(gx*r11*r21+gy*r22*r12+gz*r13*r23);
offe=-c0+wme;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hamil=zeros(2^nsp,2^nsp,nsteps);


for ii = 0:10000-1
    
    hhyp_zz=hzz_max*(-0.5*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+...
        sqrt(2)*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en))*2*IzSz;
    hhyp_zx=hzz_max*(-0.5*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en)-...
        (sqrt(2)/4)*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+(sqrt(2)/4)*(3*(cosd(beta_en)^2)-1))*(IpSz+ImSz);
   hhyp1=1*(hhyp_zz+hhyp_zx);
    ganisohamil=c0+c1*cosd(360*wr*ii*trstep)+c2*sind(360*wr*ii*trstep)+c3*cosd(360*wr*ii*trstep*2)+c4*sind(360*wr*ii*trstep*2);
    hamil(:,:,ii+1) = (ganisohamil-wme)*Sz+wn*Iz+1*hhyp1;
end

Zp=exp((we+wn)*planck/(kb*T))+exp((-we+wn)*planck/(kb*T))+exp((we-wn)*planck/(kb*T))+exp(-(we+wn)*planck/(kb*T));
rho_zeeman=(1/Zp)*expm(-(we*Sz+wn*Iz)*planck/(kb*T)); %or define it this way
ps0=trace(rho_zeeman*Sz); 
pi0=trace(rho_zeeman*Iz);
p0=pi0+ps0;  %%% Total polarisation for plotting !!!

[evecham,evalgham]=eigenshuffle(hamil); %%% Sort the eigenvalues

prop=zeros(16,16,10000);
Dinv_store=zeros(4,4,10000);
for ii=1:10000
    D=squeeze(evecham(:,:,ii));
    Dinv=D\eye(4); %%%% same as inv(D)
    hamil_mwt=Dinv*h_mw*D; %%% Transform the MW Hamiltonian to the same frame
    hamil_0=diag(evalgham(:,ii));
    hamilt=hamil_0+hamil_mwt; %%%% Total Hamiltonian
 
 %%%%%% Tilt all the operators to the same frame %%%%

    Szt=Dinv*Sz*D;
    Izt=Dinv*Iz*D;
    Sxt=Dinv*Sx*D;
    Ixt=Dinv*Ix*D;
    Syt=Dinv*Sy*D;
    Iyt=Dinv*Iy*D;
    Spt=Dinv*Sp*D;
    Smt=Dinv*Sm*D;
    Ipt=Dinv*Ip*D;
    Imt=Dinv*Im*D;
    
    
   
   %%%%%% Make left and right Liouvillian %%%%%%%%% Generate the
   %%%%%% Lindbladians
    LSzt=kron(Szt,Szt.');
    RSzt=kron(eye(2^nsp),eye(2^nsp));
    LIzt=kron(Izt,Izt.');
    RIzt=kron(eye(2^nsp),eye(2^nsp));
    LIpt=1.0*kron(Ipt,Imt.')-.5*eye(4^nsp)+.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    LImt=1.0*kron(Imt,Ipt.')-.5*eye(4^nsp)-.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    LSpt=1.0*kron(Spt,Smt.')-.5*eye(4^nsp)+.5*(kron(Szt,eye(2^nsp))+kron(eye(2^nsp),Szt.'));
    LSmt=1.0*kron(Smt,Spt.')-.5*eye(4^nsp)-.5*(kron(Szt,eye(2^nsp))+kron(eye(2^nsp),Szt.'));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Rt2ematv3=(1/t2e)*(LSzt-.25*RSzt);
    Rt2nmatv3=(1/(t2n))*(LIzt-.25*RIzt);
    Rt1mat=gep*LSpt+gem*LSmt+...%+gep*LS2pt+gem*LS2mt+...               %(-2*LSyt*RSyt+LSyt*LSyt+RSyt*RSyt)+(-2*LSzt*RSzt+LSzt*LSzt+RSzt*RSzt))+...
        gnp*LIpt+gnm*LImt;
    Rtot=1*Rt1mat+1*Rt2ematv3+1*Rt2nmatv3;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Lhamilt=kron((hamilt),eye(4))-kron(eye(4),(hamilt).');
 
    LD=kron(D,D);%-kron(eye(4),D.'); This is the transformation of propagators to the Zeeman frame; see Mance paper
    LDinv=kron(Dinv,Dinv); % or early SABRE papers Quantitative description of the SABRE process: rigorous consideration
                              %%%%of spin dynamics and chemical exchange; Konstantin Ivanov
    Ltot=Lhamilt+1*1i*((Rtot));
    prop(:,:,ii)=LD*expm(-1i*Ltot*trstep)*LDinv;
end


prop_accu=eye(16);
for kk=1:nsteps
    prop_accu=prop_accu*squeeze(prop(:,:,kk));
end


nrot=round(final_time/tr);
iznew=zeros(1,nrot);
sznew=iznew;
rho0t=rho_zeeman;


for jj=1:nrot
        iznew(jj)=trace(rho0t*Iz);
        sznew(jj)=trace(rho0t*Sz);
        Lrho0t=reshape(rho0t,[16,1]);
        rho0t=prop_accu*Lrho0t;
        rho0t=reshape(rho0t,[4,4]);
end

%toc

end
