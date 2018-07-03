clear variables
%%% Parameters *(move into dynamics file)
t1_elec = 0.3E-3;                       % Electron spin-lattice relaxtion T1 (s)
t1_nuc = 10;                            % Nuclear spin-lattice relaxation T1 (s)
freq_elec = 264e9; %265.2E9;                    % Electron Larmor frequency (Hz)
gtensor_tempol = [253.6, 105.1, 123.8]; % Tempol g-tensors
freq_rotor = 3E3;                       % Spinning frequency (Hz)
period_rotor = 1 / freq_rotor;          % Rotor period (s)
%freq_microwave = (1:5:200);             % Microwave frequency array (GHz)
freq_microwave = 0.85; %16;

final_time=40; % not sure what significance of this is
data_points = round(final_time * freq_rotor); % *Number of data points? 

%%% Pre-allocaiton
pol_iz = zeros(1, data_points); % Nuclear polorisation
pol_sz = zeros(1, data_points); % *obselete?
pol_iz_avg = zeros(length(freq_microwave), data_points); % *obselete, simply stores pol_iz?
pol_sz_avg = zeros(length(freq_microwave), data_points); % *obselete?
max_pol_iz = zeros(1, length(freq_microwave)); % Maximim nuclear polorisation

count=0;
tic
for jj=1:length(freq_microwave)
    for ii=1:length(gtensor_tempol(1)) % Powder averaging?
        
        [pol_iz(ii,:), pol_sz(ii,:), iz_rot, sz_rot, evalgham]= ... % *Why return density matrix?
            dynamics(gtensor_tempol(1), gtensor_tempol(2), gtensor_tempol(3), ...
            freq_rotor/1E3, final_time, freq_microwave(jj), ...
            freq_elec, t1_elec, t1_nuc);
        
        toc
        
        count=count+1;
        disp(['Simulation ', num2str(count), ' of ', num2str(length(freq_microwave))])
    end
    %pol_iz_avg(jj,:)=mean(pol_iz,1);
    %pol_sz_avg(jj,:)=mean(pol_sz,1);
end

toc

% for pp=1:length(freq_microwave) % *Merge loops
%     max_pol_iz(pp)=pol_iz_avg(pp, data_points);
% end
% 
% save(sprintf('%s_%s_%d.mat','solid_effect','microwave_power', ... 
%     freq_microwave(length(freq_microwave))));
% disp('End of Simulation and Data Saved')
