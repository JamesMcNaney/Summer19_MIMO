function csi_mat = channel_sim(par)

%% Set up input parameters
% feel free to change these parameters
show = 0; % 1 = generate plots, 0 = don't generate plots
% rng(1) % set random seed: comment out this line to generate a different channel each time
% par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
% par.fc = 60e9; % carrier frequency [Hz]
% par.BW = 14e6; % bandwidth [Hz]
% par.N = 2048; % number of carriers
% par.B = 64; % number of antennas in the BS (we use a single BS)
% par.U = 8; % number of single-antenna UEs

%% Set up simulation parameters
% probably don't need to change these parameters
s = qd_simulation_parameters;
s.sample_density = 2; % 2 samples per half-wavelength (spacial sampling)
s.center_frequency = par.fc;
s.use_absolute_delays = false; % Include delay of the LOS path
s.use_spherical_waves = true; % No need to use spherical waves if rx and tx are far enough
s.show_progress_bars = false;

%% Define array geometry
% feel free to change UE_x_locs, UE_y_locs, UE_z_locs, BS_x_locs, BS_y_locs, BS_z_locs

% BS is tx, UE is rx array

% UE_x_locs, UE_y_locs, UE_z_locs are column vectors We randomly place UEs
% in a rectangular  area (units in m)
% UE_x_locs = 50 + 50*rand(par.U,1);
% UE_y_locs = 50*(2*(rand(par.U,1)-0.5));
% UE_z_locs = 1.5*ones(par.U,1);
UE_x_locs = 20 + 20*rand(par.U,1);
UE_y_locs = 20*(2*(rand(par.U,1)-0.5));
UE_z_locs = 1.5*ones(par.U,1);
% place BS antennas only on y-axis at half wavelength spacing (units in m)
BS_x_locs = zeros(par.B,1);
BS_y_locs = s.wavelength/2*(-(par.B-1)/2:1:(par.B-1)/2)';
BS_z_locs = 25*ones(par.B,1); % note that the height of the transmitter is 25m

%% Assign geometry to layout object
% Create new QuaDRiGa layout object
l = qd_layout(s);

% Trajectories have only one snapshot -> UEs are static
l.track.no_snapshots = 1;
l.rx_array = qd_arrayant('omni'); % Omnidirectional UE antennas
% UEs geometry (UE is rx array)
l.no_rx = par.U; % Assign the number of UEs
l.rx_position = [UE_x_locs'; UE_y_locs'; UE_z_locs']; % Write user locations

% Define BS antenna array
l.tx_array.generate('omni'); % so far, no need to use other antenna
l.tx_array.no_elements = par.B;
% place BS antennas only on y-axis at half wavelength spacing, units in m
l.tx_array.element_position = [BS_x_locs'; BS_y_locs'; BS_z_locs'];
l.tx_position = [0;0;0]; % make l.tx_array.element_position relative to the origin

% plot scenario
if(show)
    l.visualize([],[],0); % Plot the layout
    view(-33, 60); % Enable 3D view
end

%% Set the scenario
l.set_scenario(par.scenario);
%% Generate channel coefficients
c = l.get_channels; % Generate channels
csi_mat = zeros(par.B,par.U);
for i = 1:par.U
     coll = permute(c(i,1).coeff,[2,1,3]);
     coll = sum(coll,3);
     csi_mat(:,i) = coll;
end
% csi_mat = normalize(csi_mat,1);
%%
% %% Calculating the frequency response for this channel
% 
% % Writing the frequency response in the desired way: a matrix corresponding
% % to the user grid, each point with a vector with the number of antennas
% HN = zeros(par.B,par.U,par.N);
% for uu=1:par.U
%     HN(:,uu,:) = c(uu).fr(par.BW,par.N); % Using built-in method to evaluate the frequency response from the channel coefficients
% end
% 
% if(show)
% figure(10)
% imagesc(abs(squeeze(fft(HN(:,1,:))))); colorbar
% xlabel('subcarrier')
% ylabel('spatial FFT over BS')
% 
% figure(11)
% imagesc((abs(fft(HN(:,:,1))))); colorbar
% xlabel('UE')
% ylabel('spatial FFT over BS')
% end
