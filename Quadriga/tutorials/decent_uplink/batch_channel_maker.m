% rng(21) % set random seed: comment out this line to generate a different channel each time

par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
par.fc = 3.5e9; % carrier frequency [Hz]
par.BW = 10e6; % bandwidth [Hz]
par.N = 1024; % number of carriers
par.B = 128; % number of antennas in the BS (we use a single BS)
par.U = 16; % number of single-antenna UEs
par.array_v = 4;
par.array_h = par.B/par.array_v;
par.trials = 4000;
csi_batch = zeros(par.B,par.U,par.trials);

%% Set up simulation parameters
% probably don't need to change these parameters
s = qd_simulation_parameters;
s.sample_density = 2; % 2 samples per half-wavelength (spacial sampling)
s.center_frequency = par.fc;
s.use_absolute_delays = false; % Include delay of the LOS path
s.use_spherical_waves = false; % No need to use spherical waves if rx and tx are far enough
s.show_progress_bars = false;

%%
%creating a randomization of UE's that guarantees outside of sep_ang
%degrees separation between UEs. Up to a maximum of rand_trials iterations

sep_ang = 2;                                %miminum degrees separation desired
                               %compared against rand_trials
angle_par = 1;                              %stays 1 unless UEs are adequately spaced

%%
% place BS antennas as a 2D array, same for all trials
BS_x_locs = zeros(par.array_v,par.array_h);
BS_y_locs = s.wavelength/2*(-(par.array_h-1)/2:1:(par.array_h-1)/2);
BS_z_locs = (s.wavelength/2*(-(par.array_v-1)/2:1:(par.array_v-1)/2))';

BS_y_locs = kron(ones(par.array_v,1),BS_y_locs);
BS_z_locs = 25*ones(par.array_v, par.array_h) + kron(ones(1,par.array_h),BS_z_locs);

BS_x_locs = reshape(BS_x_locs',par.B,1);
BS_y_locs = reshape(BS_y_locs',par.B,1);
BS_z_locs = reshape(BS_z_locs',par.B,1);

%%
test_graph = zeros(3,par.U,par.trials);
%%
% Create new QuaDRiGa layout object
l = qd_layout(s);

% Trajectories have only one snapshot -> UEs are static
l.track.no_snapshots = 1;

% UEs geometry (UE is rx array)
l.no_rx = par.U; % Assign the number of UEs
% Define BS antenna array
l.tx_array.generate('omni'); % so far, no need to use other antenna
l.tx_array.no_elements = par.B;
% place BS antennas only on y-axis at half wavelength spacing, units in m
l.tx_array.element_position = [BS_x_locs'; BS_y_locs'; BS_z_locs'];
l.tx_position = [0;0;0]; % make l.tx_array.element_position relative to the origin

%% Set the scenario
l.set_scenario(par.scenario);

for trial = 1:par.trials
    rng(trial);
    complete = 100*trial/par.trials;
    disp([num2str(complete,'%.3f') '% completed'])
    
UE_dist = 50 + 50*rand(par.U,1);            %distance from BS between 50 and 100
UE_ang = zeros(1,16);
while(angle_par)   %while less than rand_trials run and UEs don't pass spacing
    UE_ang = -60+120*rand(1,par.U);                                 %generate random angles (-60,60) for each UE
    check = toeplitz([UE_ang(1) fliplr(UE_ang(2:end))], UE_ang);    %make array into circulant matrix
    check = abs(check-(UE_ang.*ones(par.U,par.U)));                 %find difference between each user and all its neighbors
    check = min(min(check(2:end,:)));                               %find minimum angles between UEs (not including with itself)
    if(check>sep_ang)
        angle_par = 0;
    end
end
angle_par = 1;
UE_ang = (pi*UE_ang/180)';              %put angles into radian
UE_x_locs = UE_dist.*cos(UE_ang);       %'polar' coordinate to rectangular
UE_y_locs = UE_dist.*sin(UE_ang);       %'polar' coordinate to rectangular
UE_z_locs = 1.5*ones(par.U,1);



%% Assign geometry to layout object
l.rx_array = qd_arrayant('omni'); % Omnidirectional UE antennas
l.rx_position = [UE_x_locs'; UE_y_locs'; UE_z_locs']; % Write user locations

test_graph(:,:,trial) = [UE_x_locs'; UE_y_locs'; UE_z_locs'];
%% Generate channel coefficients
cb = l.init_builder;
cb.plpar = [];
cb.gen_ssf_parameters;
c = (cb.get_channels)'; % Generate channels for NLOS scenarios     
% c = l.get_channels; % Generate channels for LOS scenarios
csi_mat = zeros(par.B,par.U);
for i = 1:par.U
     coll = permute(c(i,1).coeff,[2,1,3]);
     coll = sum(coll,3);
     csi_mat(:,i) = coll;
end
    norm_coef = zeros(1,par.U);                    
        for i = 1:par.U
            for j = 1:par.B
                norm_coef(i)=norm_coef(i)+norm(csi_mat(j,i)); %sum the 2-norms of each column
            end
            norm_coef(i) = norm_coef(i)/par.B;         %average the 2-norm sum
            csi_mat(:,i) = csi_mat(:,i)/norm_coef(i)*.07833;
        end
    csi_batch(:,:,trial) = csi_mat;
end