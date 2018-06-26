function[iznew,sznew,ps0,pi0,evalgham]=dynamics(a,b,g,nr,final_time,mw_pwr,wme,t1e,t1n)

%%% Spin matrices (easyspin libriary)
spins = [1/2 1/2]; % S=1/2, I=1/2 system
nsp = length(spins); % Number of spins
Sz=sop(spins,'ze');
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

%%% Parameters
wn=-400.9e6; % Nuclear frequency
we=28.025*9.4*1e9; % Electron frequency (B=9.4T)
w1=mw_pwr*1e6; % Convert microwave frequency to MHz
spinning_freq=nr; % Rotor frequency
wr = spinning_freq*1000; % Rotor frequency
tr = 1/wr; % Rotor period
nsteps = 10000; % Number of timesteps (absolute)
trstep = tr/nsteps; % Value of timestep (changes with rotor period)
tarray=0:trstep:tr; % Time array
hzz_max=3e6; % Hyperfine coupling
beta_en=b; % G-tensors
gamma_en=g; % G-tensors
alpha_g=a;  % G-tensors
beta_g=b;  % G-tensors
gamma_g=g;  % G-tensors
hbar=1.054e-34; % Reduced planck
planck=hbar*2*pi; % planck
kb=1.3806e-23; % Boltzmann
T=100; % Temperature
t2e=1e-6; % T2 electron
t2n=1e-3; % T2 nucleus
pl=planck; % planck
h_mw=w1*Sx; % Microwave Hamiltonian
nrot=round(final_time/tr); % Number of data points?
iznew=zeros(1,nrot); % polarisation Iz
sznew=iznew; % polarisation Sz

% What are these, something to do with Liouvillian?
t_corr=pl/(kb*T);
p_e=tanh(0.5*we*t_corr);
p_n=tanh(0.5*wn*t_corr);
gnp=0.5*(1-p_n)*(1/(1*t1n)); % This is thermalisation, expansion of exp-(t_corr*H)/Trace(same);
gnm=0.5*(1+p_n)*(1/(1*t1n));
gep=0.5*(1-p_e)*(1/(1*t1e));
gem=0.5*(1+p_e)*(1/(1*t1e));

%%% G-anisotropy
gx=we*(2.00614/2)+18.76e6;
gy=we*(2.00194/2)+92.4e6;
gz=we*(2.00988/2)+18.2e6;

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

%%% Hamiltonian
hamil=zeros(2^nsp,2^nsp,nsteps);

for ii = 0:10000-1 % nsteps = 10000, starting from zero
    hhyp_zz=hzz_max*(-0.5*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+...
        sqrt(2)*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en))*2*IzSz;
    
%     hhyp_zz = hzz_max * (-0.5*(sind(beta_en)^2)) * IzSz;
    
    hhyp_zx=hzz_max*(-0.5*sind(beta_en)*cosd(beta_en)*cosd(360*wr*ii*trstep+gamma_en)-...
        (sqrt(2)/4)*(sind(beta_en)^2)*cosd(2*(360*wr*ii*trstep+gamma_en))+(sqrt(2)/4)*(3*(cosd(beta_en)^2)-1))*(IpSz+ImSz);
    
%     hhyp_zx = hzz_max * (IpSz+ImSz);
    
    hhyp1=1*(hhyp_zz+hhyp_zx);
    ganisohamil=c0+c1*cosd(360*wr*ii*trstep)+c2*sind(360*wr*ii*trstep)+c3*cosd(360*wr*ii*trstep*2)+c4*sind(360*wr*ii*trstep*2);
    hamil(:,:,ii+1) = (ganisohamil-wme)*Sz+wn*Iz+1*hhyp1;
end

% Density matrix
Zp=exp((we+wn)*planck/(kb*T)) + exp((-we+wn)*planck/(kb*T))+exp((we-wn)*planck/(kb*T))+exp(-(we+wn)*planck/(kb*T));
% Presumably same as (1/np.sum(boltzmann_factors)) * np.diagflat(boltzmann_factors)
rho_zeeman=(1/Zp)*expm(-(we*Sz+wn*Iz)*planck/(kb*T));
rho0t=rho_zeeman;

% disp('rho0t')
% disp(rho0t)

% Below obselete?
ps0=trace(rho_zeeman*Sz);
pi0=trace(rho_zeeman*Iz);
p0=pi0+ps0;

% Calculate and sort eigenvalues
[evecham,evalgham]=eigenshuffle(hamil); % eigenvectors, eigenvalues

% disp('hamil')
% disp(size(hamil))
% disp(hamil(:, :, 1))

% disp('evecham')
% disp(size(evecham))
% disp(evecham(:, :, 1))
% 
% disp('evalgham')
% disp(size(evalgham))
% disp(evalgham(:,1))

%%% Propogation
prop=zeros(16,16,10000);

for ii=1:10000 % nsteps = 10000
%for ii=1:1 % nsteps = 10000
    D=squeeze(evecham(:,:,ii));
    Dinv=inv(D);
    hamil_mwt=Dinv*h_mw*D; % Transform the MW Hamiltonian to the same frame
%     disp('hamil_mwt')
%     disp(hamil_mwt)
    
    hamil_0=diag(evalgham(:,ii)); % Intrinsic Hamiltonian
    hamilt=hamil_0+hamil_mwt; % Total Hamiltonian (intrinsic + MW)
%     disp('hamilt')
%     disp(hamilt)
    
    % Tilt all the operators to the same frame
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
%     
%     disp('Szt')
%     disp(Szt)
%     
%     disp('Im')
%     disp(Im)
%     
%     disp('Imt')
%     disp(Imt)
%     
%     disp('Sp')
%     disp(Sp)
%     
%     disp('Spt')
%     disp(Spt)
    
    % Liouvillian
    LSzt=kron(Szt,Szt.');
    
%     disp('LSzt')
%     disp(LSzt)
  
    RSzt=kron(eye(2^nsp),eye(2^nsp));
    LIzt=kron(Izt,Izt.');
    RIzt=kron(eye(2^nsp),eye(2^nsp)); % equal to RSzt
    
    LIpt=1.0*kron(Ipt,Imt.')-.5*eye(4^nsp)+.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
%     disp('LIpt')
%     disp(LIpt)
    
    LImt=1.0*kron(Imt,Ipt.')-.5*eye(4^nsp)-.5*(kron(Izt,eye(2^nsp))+kron(eye(2^nsp),Izt.'));
    LSpt=1.0*kron(Spt,Smt.')-.5*eye(4^nsp)+.5*(kron(Szt,eye(2^nsp))+kron(eye(2^nsp),Szt.'));
    LSmt=1.0*kron(Smt,Spt.')-.5*eye(4^nsp)-.5*(kron(Szt,eye(2^nsp))+kron(eye(2^nsp),Szt.'));
    
    Rt2ematv3=(1/t2e)*(LSzt-.25*RSzt);
    Rt2nmatv3=(1/(t2n))*(LIzt-.25*RIzt);
    Rt1mat=gep*LSpt+gem*LSmt+gnp*LIpt+gnm*LImt;
    Rtot=1*Rt1mat+1*Rt2ematv3+1*Rt2nmatv3;
    
    Lhamilt=kron((hamilt),eye(4))-kron(eye(4),(hamilt).');
%     disp('Lhamilt')
%     disp(Lhamilt)
    
    LD=kron(D,D);
    LDinv=kron(Dinv,Dinv);
    Ltot=Lhamilt+1*1i*((Rtot));
    test = LD*expm(-1i*Ltot*trstep)*LDinv;
    prop(:,:,ii)=LD*expm(-1i*Ltot*trstep)*LDinv; % Transformation back into Zeeman frame?
    
%     disp('test')
%     disp(test)
end

prop_accu=eye(16); % 16x16 identity matrix
for kk=1:nsteps % Why multiply identity matrix with prop?
    prop_accu=prop_accu*squeeze(prop(:,:,kk));
end

% Calculate polarisation? We have 10000 timesteps, so what is nrot?
for jj=1:nrot
    iznew(jj)=trace(real(rho0t)*Iz);
    sznew(jj)=trace(real(rho0t)*Sz);
    Lrho0t=reshape(rho0t,[16,1]);

%     disp('sznew(jj)')
%     disp(sznew(jj))
%     
%     disp('iznew(jj)')
%     disp(iznew(jj))
%     
%         disp('disp(Lrho0t)')
%     disp(size(Lrho0t))
%     disp(Lrho0t)
    
    rho0t=prop_accu*Lrho0t;
    
%     disp('prop_accu')
%     disp(size(prop_accu))
%     disp(prop_accu(1, 1))
%     disp(prop_accu)
    
%     disp('rho0t')
%     disp(size(rho0t))
%     disp(rho0t)
    
    rho0t=reshape(rho0t,[4,4]);
%     disp('rho0t')
%     disp(size(rho0t))
%     disp(rho0t)
    
end

end
