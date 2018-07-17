%function [iznew,sznew,evalgham] = solid_effect_dnpsims(spinning_freq)
%clear
%%%%%% Define spin operators %%%%%%
tic
spins = [1/2 1/2];
nsp = length(spins);
Sz=sop(spins,'ze');
Iz=sop(spins,'ez');
Ix=sop(spins,'ex');
Iy=sop(spins,'ey');
IxSz=sop(spins,'zx');
IySz=sop(spins,'zy');
IzSz=sop(spins,'zz');
Sp=sop(spins,'+e');
Sm=sop(spins,'-e');
Ip=sop(spins,'e+');
Im=sop(spins,'e-');
Sx=sop(spins,'xe');
Sy=sop(spins,'ye');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Define the constants and magnitudes %%%%%%%
wme=264e9;
wn=-400.9e6;
we=28.025*9.4*1e9;
%we=263e9;
%we=wme-wn;
w1=0.85e6;
delw=we-wme;
spinning_freq=3;
wr = spinning_freq*1000;
tr = 1/wr;
nsteps = 10000;
trstep = tr/nsteps;
tarray=0:trstep:tr;
hzz_max=5e6;
%rho0 = Sz;
beta_en=0;
gamma_en=0;
h_mw=w1*Sx;
hbar=1.054e-34;
planckc=hbar*2*pi;
kb=1.3806e-23;
T=100;
t1e=0.3e-3;
t1n=10;
t2e=1e-6;
t2n=1e-3;
pl=planck;
Nrot=1;
b_fact=pl/(kb*T);

t_corr=pl/(kb*T);
Bf=we*t_corr;
Bfn=wn*t_corr;
p_e=tanh(.5*we*t_corr);
p_n=tanh(.5*wn*t_corr);
gnp=0.5*(1-p_n)*(1/t1n);
gnm=0.5*(1+p_n)*(1/t1n);
gep=0.5*(1-p_e)*(1/t1e);
gem=0.5*(1+p_e)*(1/t1e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% define g-anisotropy %%%%%%%%
thetam=acosd(1/sqrt(3));
alpha_g=253.6;
beta_g=105.1;
gamma_g=123.8;
gx=we*(2.00614/2);
gy=we*(2.00194/2);
gz=we*(2.00988/2);
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

for ii = 0:length(tarray)-2
    hhyp_zz=hzz_max*(-0.5*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+...
        sqrt(2)*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en))*IzSz;
    hhyp_zx=hzz_max*(-0.5*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en)-...
        (sqrt(2)/4)*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+(sqrt(2)/4)*(3*(cosd(beta_en)^2)-1))*IxSz;
    hhyp_zy=hzz_max*((sqrt(6)/4)*(sind(beta_en)^2)*sind(2*(360*wr*ii*trstep+gamma_en))-...
        (sqrt(3)/2)*sind(beta_en)*cosd(beta_en)*sind(1*(360*wr*ii*trstep+gamma_en)))*IySz;
    hhyp=2*(hhyp_zz+hhyp_zx);
    ganisohamil=c0+c1*cosd(360*wr*ii*trstep)+c2*sind(360*wr*ii*trstep)+c3*cosd(360*wr*ii*trstep*2)+c4*sind(360*wr*ii*trstep*2);
    hamil(:,:,ii+1) = (ganisohamil-wme)*Sz+wn*Iz+hhyp;
end


Zp=exp((we+wn)*planck/(kb*T))+exp((-we+wn)*planck/(kb*T))+exp((we-wn)*planck/(kb*T))+exp(-(we+wn)*planck/(kb*T));
%rho_zeeman=(1/Zp)*expm((-squeeze(hamil(:,:,1))-wme*Sz)*planck/(kb*T)); %%% Define zeeman rho0
rho_zeeman=(1/Zp)*expm(-(we*Sz+wn*Iz)*planck/(kb*T)); %or define it this way
%rho_zeeman=Sz;   %%%%% Just for checking
ps0=trace(rho_zeeman*Sz); 
pi0=trace(rho_zeeman*Iz);
p0=pi0+ps0;  %%% Total polarisation for plotting !!!

%  Rcorre=[sum((planck/(2*kb*T))*we*(kron(Szt,eye(4))-kron(eye(4),Szt.')),2) zeros(16,15)];
%      Rcorrn=[sum((planck/(2*kb*T))*wn*(kron(Izt,eye(4))-kron(eye(4),Izt.')),2) zeros(16,15)];
%      Rcorr=Rcorre-0*Rcorrn;
%      Rcorr(:,2:end)=0;
%      Rcorr(1,1:end)=0;
% 

[evecham,evalgham]=eigenshuffle(hamil); %%% Sort the eigenvalues
% iznew=zeros(1,length(tarray)); 
% sznew=iznew;
iznew_accu=zeros(1,Nrot*length(tarray));
sznew_accu=iznew_accu;
D=zeros(4,4);  
Dinv=D;

for ii=1:length(tarray)-1
 
    D=squeeze(evecham(:,:,ii));
    Dinv=D\eye(4); %%%% same as inv(D)
    hamil_mwt=Dinv*h_mw*D; %%% Transform the MW Hamiltonian to the same frame
    hamil_0=diag(evalgham(:,ii)); 
    hamilt=hamil_0+hamil_mwt; %%%% Total Hamiltonian
   %rho_zeeman=(expm(-hamil_0)*planck/kb*T);
    eigval=evalgham(:,ii);
   %%%%%% Tilt all the operators to the same frame %%%%   
    Szt=Dinv*Sz*D;
    Izt=Dinv*Iz*D;
    Sxt=Dinv*Sx*D;
    Ixt=Dinv*Ix*D;
    Syt=Dinv*Sy*D;
    Iyt=Dinv*Iy*D;
    Spt=Dinv*Sp*D;
    Ipt=Dinv*Ip*D;
    Smt=Dinv*Sm*D;
    Imt=Dinv*Im*D;
    

   %%%%%% Make left and right Liouvillian %%%%%%%%%
    LSzt=kron(Szt,eye(4));
    RSzt=kron(eye(4),Szt');
    LIzt=kron(Izt,eye(4));
    RIzt=kron(eye(4),Izt');
    LSxt=kron(Sxt,eye(4));
    RSxt=kron(eye(4),Sxt');
    LIxt=kron(Ixt,eye(4));
    RIxt=kron(eye(4),Ixt');
    LIpt=kron(Ipt,eye(4));
    RIpt=kron(eye(4),Ipt.');
    LSpt=kron(Spt,eye(4));
    RSpt=kron(eye(4),Spt.');
    LIyt=kron(Iyt,eye(4));
    RIyt=kron(eye(4),Iyt.');
    LSyt=kron(Syt,eye(4));
    RSyt=kron(eye(4),Syt.');
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   
%%%%% Create the T1 Liouvillian %%%%%%%%%%%   
R1=zeros(2^nsp);   
R2=R1;
W2=R1;
for i=1:4
    for j=1:4
        ep=Szt(i,i)-Szt(j,j);
        epn=Izt(i,i)-Izt(j,j);
        Ef=exp(-ep*Bf/2-epn*Bfn/2)/(exp(ep*Bf/2+epn*Bfn/2)+exp(-ep*Bf/2-epn*Bfn/2));
        
        R1(i,j)=(1/t1e)*(Sxt(i,j)*Sxt(j,i)+Szt(i,j)*Szt(j,i))+ ...
                (1/t1n)*(Ixt(i,j)*Ixt(j,i)+Izt(i,j)*Izt(j,i));
            
        R1(i,j)=R1(i,j)*Ef;
    end
end

for i=1:4
    R1(i,i)=0;
end

for i=1:4
    for j=1:4
    if (abs(i-j)>0)
        R1(i,i)=R1(i,i)-R1(j,i);
    end
    end
end
% W1=expm(R1*trstep);   

Rfull=zeros(16,16);

for ib=1:4
    for jb=1:4
   
        II=(ib-1)*4+ib;
        JJ=(jb-1)*4+jb;
        IJ=(ib-1)*4+jb;
        
    Rfull(II,JJ)=R1(ib,jb);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Create T2 Liouvillian

for mm=2:4
        t2val_v2(mm,1)=abs(Szt(mm,mm)-Szt(mm-1,mm-1))^2*(1/t2e)...
                        +abs(Izt(mm,mm)-Izt(mm-1,mm-1))^2*(1/t2n);
end
relmatt2_temp=[t2val_v2; circshift(t2val_v2,1); circshift(t2val_v2,2); circshift(t2val_v2,3)];
relmatt2=diag(relmatt2_temp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Rtot=-1*relmatt2+Rfull;
   
    Lhamilt=kron((hamilt),eye(4))-kron(eye(4),(hamilt).');
 
    LD=kron(D,D);%-kron(eye(4),D.');
    LDinv=kron(Dinv,Dinv);%-kron(eye(4),Dinv.');
    Ltot=Lhamilt+1*1i*((Rtot));
    prop(:,:,ii)=LD*expm(-1i*Ltot*trstep)*LDinv;
end



prop_accu=eye(16);
for kk=1:nsteps
    prop_accu=prop_accu*squeeze(prop(:,:,kk));
end


final_time=100;
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

rho0t=rho_zeeman;

for ii=1:nsteps
        iz_rot(ii)=trace(rho0t*Iz); 
        sz_rot(ii)=trace(rho0t*Sz);
        Lrho0t=reshape(rho0t,[16,1]);
        rho0t=squeeze(prop(:,:,ii))*Lrho0t;
        rho0t=reshape(rho0t,[4,4]); 
end

xbup=0+tr/2:tr:final_time-tr/2;
xrot=(0+trstep/2:trstep:tr-trstep/2).*1e6;

figure;
subplot(3,1,1)
plot(xrot,real(iz_rot/pi0));ylabel('P_i/P_{i0}')
subplot(3,1,2)
plot(xrot,real(sz_rot/ps0));ylabel('P_s/P_{s0}')
subplot(3,1,3)
plot(xrot,evalgham.');ylabel('Energy')
xlabel('Time(\mus)')

figure;
subplot(2,1,1)
plot(xbup,real(iznew/pi0));ylabel('P_i/P_{i0}')
subplot(2,1,2)
plot(xbup,real(sznew/ps0));ylabel('P_s/P_{s0}')
xlabel('Time(s)')

toc