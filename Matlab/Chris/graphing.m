clear variables %, close all
%load('solid_effect_microwave_power_16.mat');
load('solid_effect_microwave_power_8.500000e-01.mat');

% time = linspace(0,length(pol_iz_avg)-1,length(pol_iz_avg));
% enhancement_max = abs(max_pol_iz)/abs(max_pol_iz(1));
% enhancement = abs(pol_iz_avg);
% elec = abs(pol_sz_avg);

tr = 1/3E3;
final_time=1;
nsteps = 10000;
trstep = tr/nsteps;

xbup=0+tr/2:tr:final_time-tr/2;
xrot=(0+trstep/2:trstep:tr-trstep/2).*1e6;
xrot=linspace(0, 360, 10000);

figure;
subplot(3,1,1)
plot(xrot,real(iz_rot/iz_rot(1)));
ylabel('P_i/P_{i0}')
subplot(3,1,2)
plot(xrot,real(sz_rot/sz_rot(1)));
ylabel('P_s/P_{s0}')
subplot(3,1,3)
plot(xrot,evalgham);
ylabel('Energy')
xlabel('Time(\mus)')

% fig1=figure();
% plot(freq_microwave, enhancement_max, 'k', freq_microwave, enhancement_max, 'kx');
% xlabel('MW power')
% ylabel('Enhancement')
% 
% fig2=figure();
% hold on;
% plot(time, enhancement(1, :)/abs(max_pol_iz(1)), 'k', 'DisplayName', ...
%      num2str(freq_microwave(1)));
% plot(time, enhancement(5, :)/abs(max_pol_iz(1)), 'r', 'DisplayName',...
%      num2str(freq_microwave(5)));
% plot(time, enhancement(10, :)/abs(max_pol_iz(1)), 'g', 'DisplayName',...
%      num2str(freq_microwave(10)));
% plot(time, enhancement(20, :)/abs(max_pol_iz(1)), 'b', 'DisplayName',...
%      num2str(freq_microwave(20)));
 
%  plot(time, enhancement(1, :), 'k', 'DisplayName', ...
%      num2str(freq_microwave(1)));
% plot(time, enhancement(5, :), 'r', 'DisplayName',...
%      num2str(freq_microwave(5)));
% plot(time, enhancement(10, :), 'g', 'DisplayName',...
%      num2str(freq_microwave(10)));
% plot(time, enhancement(20, :), 'b', 'DisplayName',...
%      num2str(freq_microwave(20)));
% hold off;
% lgd = legend();
% title(lgd,'Microwave amplitude')
% xlabel('Time?')
% ylabel('Nuclear enhancement')