%%% Parameters
t1_elec = 0.3E-3;                       % Electron spin-lattice relaxtion T1 (s)
t1_nuc = 10;                            % Nuclear spin-lattice relaxation T1 (s)
freq_elec = 265.2E9;                    % Electron Larmor frequency (Hz)
gtensor_tempol = [253.6, 105.1, 123.8]; % Tempol g-tensors
freq_rotor = 3E3;                       % Spinning frequency (Hz)
period_rotor = 1 / freq_rotor;          % Rotor period (s)
freq_microwave = (1:5:20);              % Microwave frequency array (GHz)
data_points = round(length(freq_microwave) * freq_rotor); % Number of data points?

%%% Pre-allocaiton
pol_iz = zeros(1, data_points); % Nuclear polorisation
pol_sz = zeros(1, data_points); % Electron polorisation
ps = zeros(1, length(freq_microwave)); % ?? 
pi = zeros(1, length(freq_microwave)); % ??
pol_iz_avg = zeros(length(freq_microwave), data_points);
pol_sz_avg = zeros(length(freq_microwave), data_points);
max_pol_iz = zeros(1, length(freq_microwave));

count=0;
for jj=1:length(freq_microwave)
    for ii=1:length(gtensor_tempol(1)) % Why loop over g tensor alpha values?
        
        [pol_iz(ii,:), pol_sz(ii,:), ps(jj), pi(jj)]= ...
            fn_se_buildup(gtensor_tempol(1), gtensor_tempol(2), gtensor_tempol(3), ...
            freq_rotor/1E3, length(freq_microwave), freq_microwave(jj), ...
            freq_elec, t1_elec, t1_nuc);
        
        count=count+1;
        disp(['Simulation ', num2str(count), ' of ', ...
            num2str(length(gtensor_tempol(1))*length(freq_microwave))])
    end
    pol_iz_avg(jj,:)=mean(pol_iz,1);
    pol_sz_avg(jj,:)=mean(pol_sz,1);
end

for pp=1:length(freq_microwave)
    max_pol_iz(pp)=pol_iz_avg(pp, data_points);
end

save(sprintf('%s_%s.mat','solid_effect','microwave_power'));
disp('End of Simulation and Data Saved')
