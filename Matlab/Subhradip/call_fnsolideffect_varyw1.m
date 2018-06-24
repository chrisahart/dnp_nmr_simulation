tic
%[a,b,g]=genZCW(3);
a=253.6;
b=105.1;
g=123.8;
final_time=40;
%iz=zeros(length(a_g),Nrot*10000);
%sz=iz;
t1e=0.3*1e-3;
t1n=10;
nr=3;
nrot=round(final_time/(1/(nr*1000)));
% iz3=zeros(length(a),nrot);
% sz3=iz3;
mw_pwr=(1:5:200); %Please put the MW amplitude in MHz
wme=265.2e9;

% iz3_avg=zeros(length(wme),nrot);
% sz3_avg=iz3_avg;
count=0;
for jj=1:length(mw_pwr)
for ii=1:length(a)
    [iz3(ii,:),sz3(ii,:),ps03(jj),pi03(jj)]=...
        fn_se_buildup(a(ii),b(ii),g(ii),nr,final_time,mw_pwr(jj),wme,t1e,t1n);
    count=count+1;
    X=['Simulation ', num2str(count), ' of ', num2str(length(a)*length(mw_pwr))];
    disp(X)    
end
iz3_avg(jj,:)=mean(iz3,1);
sz3_avg(jj,:)=mean(sz3,1);
end

tr=1/(nr*1000);
xaxis=0+tr/2:tr:final_time-tr/2;

for pp=1:length(mw_pwr)
    max_iz(pp)=iz3_avg(pp,round(final_time/tr));
end

type='solid_effect';
vary_param='microwave_power';
d=datestr(datetime);
savename = sprintf('%s_%s.mat',type,vary_param);
save(savename);
Y='End of Simulation and Data Saved';
disp(Y)
toc
