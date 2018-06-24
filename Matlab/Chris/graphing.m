clear all
load('solid_effect_microwave_power.mat');

time = linspace(0,length(pol_iz_avg)-1,length(pol_iz_avg));
enhancement_max = abs(max_pol_iz)/abs(max_pol_iz(1));
enhancement = abs(pol_iz_avg);

figure(1)
plot(mw_pwr, enhancement_max, 'k', mw_pwr, enhancement_max, 'kx');
xlabel('MW power')
ylabel('Enhancement')

figure(2)
hold on;
plot(time, enhancement(5, :)/enhancement(5, 1), 'r', 'DisplayName','5');
plot(time, enhancement(10, :)/enhancement(10, 1), 'g', 'DisplayName','10');
plot(time, enhancement(20, :)/enhancement(20, 1), 'b', 'DisplayName','20');
hold off;
legend()
xlabel('Time?')
ylabel('Enhancement')