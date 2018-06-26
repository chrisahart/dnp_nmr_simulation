clear variables, close all
load('solid_effect_microwave_power_16.mat');

time = linspace(0,length(pol_iz_avg)-1,length(pol_iz_avg));
enhancement_max = abs(max_pol_iz)/abs(max_pol_iz(1));
enhancement = abs(pol_iz_avg);
elec = abs(pol_sz_avg);

% fig1=figure();
% plot(freq_microwave, enhancement_max, 'k', freq_microwave, enhancement_max, 'kx');
% xlabel('MW power')
% ylabel('Enhancement')

fig2=figure();
plot(abs(pol_iz(1, :)));

fig3=figure();
plot(evalgham(1, :), 'r')
hold on;
plot(evalgham(2, :), 'g')
plot(evalgham(3, :), 'b')
plot(evalgham(4, :), 'y')
hold off;

% fig2=figure();
% hold on;
% plot(time, enhancement(5, :)/enhancement(5, 1), 'r', 'DisplayName','5');
% plot(time, enhancement(10, :)/enhancement(10, 1), 'g', 'DisplayName','10');
% plot(time, enhancement(20, :)/enhancement(20, 1), 'b', 'DisplayName','20');
% hold off;
% legend()
% xlabel('Time?')
% ylabel('Nuclear enhancement')
% 
% elec_end = 10;
% fig3=figure();
% hold on;
% plot(time(1:elec_end), elec(5, 1:elec_end)/elec(5, 1), 'r', 'DisplayName','5');
% plot(time(1:elec_end), elec(10, 1:elec_end)/elec(10, 1), 'g', 'DisplayName','10');
% plot(time(1:elec_end), elec(20, 1:elec_end)/elec(20, 1), 'b', 'DisplayName','20');
% hold off;
% legend()
% xlabel('Time?')
% ylabel('Electron polarisation')