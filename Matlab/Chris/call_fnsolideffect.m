function call_fnsolideffect()
tic
[a,b,g]=genZCW(10);
%a=253.6;
%b=105.1;
%g=123.8;
final_time=200;
%iz=zeros(length(a_g),Nrot*10000);
%sz=iz;
t1e=0.3*1e-3;
t1n=[1 4 10];
nr=3;
nrot=round(final_time/(1/(nr*1000)));
iz=zeros(length(a),nrot);
sz=iz;
mw_pwr=0.85; %Please put the MW amplitude in MHz
wme=264*1e9;

iz_avg=zeros(length(t1n),nrot);
sz_avg=iz_avg;
ps0=zeros(1,length(t1n));
pi0=ps0;
for jj=1:length(t1n)
for ii=1:length(a)
    [iz(ii,:),sz(ii,:),ps0(jj),pi0(jj)]=...
        fn_se_buildup(a(ii),b(ii),g(ii),nr,final_time,mw_pwr,wme,t1e,t1n(jj));
end
iz_avg(jj,:)=mean(iz,1);
sz_avg(jj,:)=mean(sz,1);
end

tr=1/(nr*1000);
xaxis=0+tr/2:tr:final_time-tr/2;

% subplot(1,2,1)
% plot(xaxis(1:1000:end),abs(iz_avg(1:1000:end)/pi0));
% subplot(1,2,2)
% plot(xaxis(1:1000:end),abs(iz_avg(1:1000:end)/ps0));

save solid_effect_bup_17Jun2018_varyt1n_t1econst_timeindhyp.mat
toc
end
