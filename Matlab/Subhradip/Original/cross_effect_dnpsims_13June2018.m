%function [iznew,sznew,evalgham] = solid_effect_dnpsims(spinning_freq)
tic
clear
%%%%%% Define spin operators %%%%%%
spins = [1/2 1/2 1/2];
nsp = length(spins);
S1z=sop(spins,'zee');
S2z=sop(spins,'eze');
Iz=sop(spins,'eez');
Ix=sop(spins,'eex');
Iy=sop(spins,'eey');
S1zIx=sop(spins,'zex');
S1zIy=sop(spins,'zey');
S1zIz=sop(spins,'zez');
S2zIx=sop(spins,'ezx');
S2zIy=sop(spins,'ezy');
S2zIz=sop(spins,'ezz');
S1x=sop(spins,'xee');
S1y=sop(spins,'yee');
S2x=sop(spins,'exe');
S2y=sop(spins,'eye');
S1zS2z=sop(spins,'zze');
S1pS2m=sop(spins,'+-e');
S1mS2p=sop(spins,'-+e');
S20=2*S1zS2z-1/2*(S1pS2m+S1mS2p);
Ip=sop(spins,'ee+');
Im=sop(spins,'ee-');
S1p=sop(spins,'+ee');
S1m=sop(spins,'-ee');
S2p=sop(spins,'e+e');
S2m=sop(spins,'e-e');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Define the constants and magnitudes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr=3;
final_time=100;
wme=264e9;
wn=-400.9e6;
we=28.025*9.4*1e9;
mw_pwr=0.85;
w1=mw_pwr*1e6;
spinning_freq=nr;
wr = spinning_freq*1000;
tr = 1/wr;
nrot=round(final_time/tr);
nsteps = 10000;
trstep = tr/nsteps;
hzz_max=3e6;
beta_en=0;
gamma_en=0;
h_mw=w1*(S1x+S2x);
kb=1.3806e-23;
T=100;
t1n=10;
t1e=0.3e-3;
t2e=1e-6;
t2n=1e-3;
pl=planck;
%t_corr=pl/(kb*T);
t_corr=pl/(kb*T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% define g-anisotropy %%%%%%%%
alpha_g=253.6;
beta_g=105.1;
gamma_g=123.8;
[c01,c11,c21,c31,c41]=ganiso(alpha_g,beta_g,gamma_g,we);
[c02,c12,c22,c32,c42]=ganiso(alpha_g+102,beta_g+104,gamma_g+124,we);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hamil=zeros(2^nsp,2^nsp,nsteps);
b_ee=0;
g_ee=0;
for ii = 0:nsteps-1
    hhyp_zz=hzz_max*(-0.5*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+...
        sqrt(2)*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en))*S1zIz;
    hhyp_zx=hzz_max*(-0.5*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en)-...
        (sqrt(2)/4)*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+(sqrt(2)/4)*(3*(cosd(beta_en)^2)-1))*S1zIx;
    hhyp_zy=hzz_max*((sqrt(6)/4)*(sind(beta_en)^2)*sind(2*(360*wr*ii*trstep+gamma_en))-...
        (sqrt(3)/2)*sind(beta_en)*cosd(beta_en)*sind(1*(360*wr*ii*trstep+gamma_en)))*S1zIy;
    hhyp=2*(hhyp_zz+hhyp_zx+0*hhyp_zy);
    ganisohamil1=c01+c11*cosd(360*wr*ii*trstep)+c21*sind(360*wr*ii*trstep)+c31*cosd(360*wr*ii*trstep*2)...
        +c41*sind(360*wr*ii*trstep*2);
    ganisohamil2=c02+c12*cosd(360*wr*ii*trstep)+c22*sind(360*wr*ii*trstep)+c32*cosd(360*wr*ii*trstep*2)...
        +c42*sind(360*wr*ii*trstep*2);
    dipe=23e6*(+0.5*((sind(b_ee))^2)*cosd(2*360*wr*ii*trstep+g_ee)-...
        sqrt(2)*sind(b_ee)*cosd(b_ee)*cosd(360*wr*ii*trstep+g_ee))*S20;
    hamil(:,:,ii+1) = (ganisohamil1-wme)*S1z+(ganisohamil2-wme)*S2z+wn*Iz+1*hhyp+0*dipe;
end


elev_h0=diag(squeeze(hamil(:,:,1)));
Zp=sum(exp(-(elev_h0*planck)/(kb*T)));
rho_zeeman=(1/Zp)*expm(-(we*S1z+we*S2z+1*wn*Iz)*planck/(kb*T)); %or define it this way
%rho_zeeman=S1z+S2z;
ps10=trace(rho_zeeman*S1z);
ps20=trace(rho_zeeman*S2z); 
pi0=trace(rho_zeeman*Iz);
p0=pi0+ps10+ps20;  
[evecham,evalgham]=eigenshuffle(hamil); %%% Sort the eigenvalues
gnp=1/t1n*exp(-0.5*t_corr*wn)/Zp;
gnm=1/t1n*exp(0.5*t_corr*wn)/Zp;
gep=1/t1e*exp(-0.5*t_corr*we)/Zp;
gem=1/t1e*exp(0.5*t_corr*we)/Zp;


prop=zeros(64,64,nsteps);
for ii=1:nsteps
    D=squeeze(evecham(:,:,ii));
    Dinv=D\eye(8); %%%% same as inv(D)
    hamil_mwt=Dinv*h_mw*D; %%% Transform the MW Hamiltonian to the same frame
    hamil_0=diag(evalgham(:,ii));
    hamilt=hamil_0+hamil_mwt; %%%% Total Hamiltonian
 
 %%%%%% Tilt all the operators to the same frame %%%%

    S1zt=Dinv*S1z*D;
    Izt=Dinv*Iz*D;
    S1xt=Dinv*S1x*D;
    Ixt=Dinv*Ix*D;
    S1yt=Dinv*S1y*D;
    Iyt=Dinv*Iy*D;
    S2zt=Dinv*S2z*D;
    S2xt=Dinv*S2x*D;
    S2yt=Dinv*S2y*D;
    S1pt=Dinv*S1p*D;
    S2pt=Dinv*S2p*D;
    S1mt=Dinv*S1m*D;
    S2mt=Dinv*S2m*D;
    Ipt=Dinv*Ip*D;
    Imt=Dinv*Im*D;
    
    
%     S1zt=S1z;
%     Izt=Iz;
%     S1xt=S1x;
%     Ixt=Ix;
%     S1yt=S1y;
%     Iyt=Iy;
%     S2zt=S2z;
%     S2xt=S2x;
%     S2yt=S2y;
%     S1pt=S1p;
%     S2pt=S2p;
%     S1mt=S1m;
%     S2mt=S2m;
%     Ipt=Ip;
%     Imt=Im;
%     
   
   %%%%%% Make left and right Liouvillian %%%%%%%%% Generate the
   %%%%%% Lindbladians
    LS1zt=kron(S1zt,S1zt.');
    RS1zt=kron(eye(2^nsp),eye(2^nsp));
    LS2zt=kron(S2zt,S2zt.');
    RS2zt=kron(eye(2^nsp),eye(2^nsp));
    LIzt=kron(Izt,Izt.');
    RIzt=kron(eye(2^nsp),eye(2^nsp));
    
    LIpt=1.0*kron(Ipt,Imt.')-.5*eye(4^nsp)+.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    LImt=1.0*kron(Imt,Ipt.')-.5*eye(4^nsp)-.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    
    LS1pt=1.0*kron(S1pt,S1mt.')-.5*eye(4^nsp)+.5*(kron(S1zt,eye(2^nsp))+kron(eye(2^nsp),S1zt.'));
    LS1mt=1.0*kron(S1mt,S1pt.')-.5*eye(4^nsp)-.5*(kron(S1zt,eye(2^nsp))+kron(eye(2^nsp),S1zt.'));
    
    LS2pt=1.0*kron(S2pt,S2mt.')-.5*eye(4^nsp)+.5*(kron(S2zt,eye(2^nsp))+kron(eye(2^nsp),S2zt.'));
    LS2mt=1.0*kron(S2mt,S2pt.')-.5*eye(4^nsp)-.5*(kron(S2zt,eye(2^nsp))+kron(eye(2^nsp),S2zt.'));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Rt2ematv3=(2/t2e)*(LS1zt-.25*RS1zt);
    Rt2e2matv3=(2/t2e)*(LS2zt-.25*RS2zt);   
    Rt2nmatv3=(2/(t2n))*(LIzt-.25*RIzt);
    Rt1mat=gep*LS1pt+gem*LS1mt+gep*LS2pt+gem*LS2mt+...               %(-2*LSyt*RSyt+LSyt*LSyt+RSyt*RSyt)+(-2*LSzt*RSzt+LSzt*LSzt+RSzt*RSzt))+...
        gnp*LIpt+gnm*LImt;
    Rtot=1*Rt1mat+1*Rt2ematv3+1*Rt2e2matv3+1*Rt2nmatv3;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     Rt2ematv3=(1/t2e)*(LS1zt*LS1zt-2*LS1zt*RS1zt+RS1zt*RS1zt);
%     Rt2e2matv3=(1/t2e)*(LS2zt*LS2zt-2*LS2zt*RS2zt+RS2zt*RS2zt);   
%     Rt2nmatv3=(1/t2n)*(LIzt*LIzt-2*RIzt*LIzt+RIzt*RIzt);
%     Rt2ematv3=(1/t2e)*(LS1zt-0.25*RS1zt);
%     Rt2e2matv3=(1/t2e)*(LS2zt-0.25*RS2zt);   
%     Rt2nmatv3=(1/t2n)*(LIzt-0.25*RIzt);
%     Rt1mat=gep*LS1pt+gem*LS1mt+gep*LS2pt+gem*LS2mt+...               %(-2*LSyt*RSyt+LSyt*LSyt+RSyt*RSyt)+(-2*LSzt*RSzt+LSzt*LSzt+RSzt*RSzt))+...
%         gnp*LIpt+gnm*LImt;
%     Rtot=1*Rt1mat+1*Rt2ematv3+1*Rt2e2matv3+1*Rt2nmatv3;

    Lhamilt=kron((hamilt),eye(8))-kron(eye(8),(hamilt).');
 
    LD=kron(D,D);%-kron(eye(4),D.');
    LDinv=kron(Dinv,Dinv);%-kron(eye(4),Dinv.');
    Ltot=Lhamilt+1*1i*((Rtot));
    prop(:,:,ii)=LD*expm(-1i*Ltot*trstep)*LDinv;
    
%     t=expm(-1i*Ltot*trstep);
%     tic
%     for test=1:1E4
%         %t=expm(-1i*Ltot*trstep);
%         %LD=kron(D,D);
%         %Lhamilt=kron((hamilt),eye(8));
%         %S1zt=Dinv*S1z*D;
%         prop(:,:,ii)=LD*t;
%     end
%     toc
    
end

prop_accu=eye(64);
for kk=1:nsteps
    prop_accu=prop_accu*squeeze(prop(:,:,kk));
end

        

iz_rot=zeros(1,nsteps);
s1z_rot=iz_rot;
s2z_rot=iz_rot;




iz_bup=zeros(1,nrot);
s1z_bup=iz_bup;
s2z_bup=iz_bup;
rho0t=rho_zeeman;
for jj=1:nrot
        iz_bup(jj)=trace(rho0t*Iz);
        s1z_bup(jj)=trace(rho0t*S1z);
        s2z_bup(jj)=trace(rho0t*S2z);
        Lrho0t=reshape(rho0t,[64,1]);
        rho0t=prop_accu*Lrho0t;
        rho0t=reshape(rho0t,[8,8]);
        if jj==round(20/tr)
            rho0tstore=rho0t;
        end
end


rho0t=rho_zeeman;

for ii=1:nsteps
        iz_rot(ii)=trace(rho0t*Iz); 
        s1z_rot(ii)=trace(rho0t*S1z);
        s2z_rot(ii)=trace(rho0t*S2z);
        Lrho0t=reshape(rho0t,[64,1]);
        rho0t=squeeze(prop(:,:,ii))*Lrho0t;
        rho0t=reshape(rho0t,[8,8]); 
end
%end


xbup=0+tr/2:tr:final_time-tr/2;
xrot=(0+trstep/2:trstep:tr-trstep/2).*1e6;

figure;
subplot(3,1,1)
plot(xrot,real(iz_rot/-pi0));ylabel('P_i/P_{s0}')
subplot(3,1,2)
plot(xrot,real(s1z_rot/ps10));ylabel('P_s/P_{s0}')
subplot(3,1,3)
plot(xrot,evalgham.');ylabel('Energy')
xlabel('Time(\mus)')

figure;
%subplot(2,1,1)
plot(xbup,real(-iz_bup/ps10));ylabel('P_i/P_{i0}')
hold on
plot(xbup,real(s1z_bup/ps10));ylabel('P_s/P_{s0}')
plot(xbup,real(s2z_bup/ps10));ylabel('P_s/P_{s0}')
xlabel('Time(s)')

toc

