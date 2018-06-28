clear variables, close all
load('solid_effect_microwave_power_196.mat');

time = linspace(0,length(pol_iz_avg)-1,length(pol_iz_avg));
enhancement_max = abs(max_pol_iz)/abs(max_pol_iz(1));
enhancement = abs(pol_iz_avg);
elec = abs(pol_sz_avg);

fig1=figure();
plot(freq_microwave, enhancement_max, 'k', freq_microwave, enhancement_max, 'kx');
xlabel('MW power')
ylabel('Enhancement')

fig2=figure();
hold on;
% plot(time, enhancement(1, :)/abs(max_pol_iz(1)), 'k', 'DisplayName', ...
%      num2str(freq_microwave(1)));
% plot(time, enhancement(5, :)/abs(max_pol_iz(1)), 'r', 'DisplayName',...
%      num2str(freq_microwave(5)));
% plot(time, enhancement(10, :)/abs(max_pol_iz(1)), 'g', 'DisplayName',...
%      num2str(freq_microwave(10)));
% plot(time, enhancement(20, :)/abs(max_pol_iz(1)), 'b', 'DisplayName',...
%      num2str(freq_microwave(20)));
 
 plot(time, enhancement(1, :), 'k', 'DisplayName', ...
     num2str(freq_microwave(1)));
plot(time, enhancement(5, :), 'r', 'DisplayName',...
     num2str(freq_microwave(5)));
plot(time, enhancement(10, :), 'g', 'DisplayName',...
     num2str(freq_microwave(10)));
plot(time, enhancement(20, :), 'b', 'DisplayName',...
     num2str(freq_microwave(20)));
hold off;
lgd = legend();
title(lgd,'Microwave amplitude')
xlabel('Time?')
ylabel('Nuclear enhancement')