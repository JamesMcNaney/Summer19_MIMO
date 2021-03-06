function csi_mat = test04(par)

% par.array_v = 1;
% par.array_h = 128;
%% Set up input parameters
% feel free to change these parameters
show = 0; % 1 = generate plots, 0 = don't generate plots
% rng(1) % set random seed: comment out this line to generate a different channel each time
% par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
% par.fc = 3.5e9; % carrier frequency [Hz]
% par.BW = 10e6; % bandwidth [Hz]
% par.N = 1024; % number of carriers
% par.B = 128; % number of antennas in the BS (we use a single BS)
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
% UE_x_locs = 20 + 20*rand(par.U,1);
% UE_y_locs = 20*(2*(rand(par.U,1)-0.5));
% UE_z_locs = 1.5*ones(par.U,1);

%%
%creating a randomization of UE's that guarantees outside of sep_ang
%degrees separation between UEs. Up to a maximum of rand_trials iterations

sep_ang = 5;                                %miminum degrees separation desired
rand_trials = 100;                          %number of attempts to randomly create angular spacing
count = 0;                                  %compared against rand_trials
angle_par = 1;                              %stays 1 unless UEs are adequately spaced

UE_dist = 50 + 50*rand(par.U,1);            %distance from BS between 50 and 100

while((count < rand_trials) && angle_par)   %while less than rand_trials run and UEs don't pass spacing
    UE_ang = -60+120*rand(1,par.U);                                 %generate random angles (-60,60) for each UE
    check = toeplitz([UE_ang(1) fliplr(UE_ang(2:end))], UE_ang);    %make array into circulant matrix
    check = abs(check-(UE_ang.*ones(par.U,par.U)));                 %find difference between each user and all its neighbors
    check = min(min(check(2:end,:)));                               %find minimum angles between UEs (not including with itself)
    if(check>sep_ang)
        angle_par = 0;
    else
        count = count+1;
    end
end
%If the randomly generated angles fail to meet the sep_ang criteria, print
%out to the command window, but still carry on with program.
if count == rand_trials
    'Angle miminum between UEs NOT met'             
end
UE_ang = (pi*UE_ang/180)';              %put angles into radian
UE_x_locs = UE_dist.*cos(UE_ang);       %'polar' coordinate to rectangular
UE_y_locs = UE_dist.*sin(UE_ang);       %'polar' coordinate to rectangular
z_r = randi([0,5],par.U,1);
UE_z_locs = 4*z_r + 1.5;

% place BS antennas only on y-axis at half wavelength spacing (units in m)
% BS_x_locs = zeros(par.B,1);
% BS_y_locs = s.wavelength/2*(-(par.B-1)/2:1:(par.B-1)/2)';
% BS_z_locs = 25*ones(par.B,1); % note that the height of the transmitter is 25m

% place BS antennas as a 2D array
BS_x_locs = zeros(par.array_v,par.array_h);
BS_y_locs = s.wavelength/2*(-(par.array_h-1)/2:1:(par.array_h-1)/2);
BS_z_locs = (s.wavelength/2*(-(par.array_v-1)/2:1:(par.array_v-1)/2))';

BS_y_locs = kron(ones(par.array_v,1),BS_y_locs);
BS_z_locs = 25*ones(par.array_v, par.array_h) + kron(ones(1,par.array_h),BS_z_locs);

BS_x_locs = reshape(BS_x_locs',par.B,1);
BS_y_locs = reshape(BS_y_locs',par.B,1);
BS_z_locs = reshape(BS_z_locs',par.B,1);


%% Assign geometry to layout object
% Create new QuaDRiGa layout object
l = qd_layout(s);

% Trajectories have only one snapshot -> UEs are static
% l.track.no_snapshots = 1;

l.track(1,1) = qd_track('circular',20*pi);
l.track(1,1).name = 'Rx1';                              % Set the MT1 name

l.track(1,2) = qd_track('circular',20*pi);
l.track(1,2).name = 'Rx2';                              % Set the MT1 name

l.track(1,3) = qd_track('circular',20*pi);
l.track(1,3).name = 'Rx3';                              % Set the MT1 name

l.track(1,4) = qd_track('circular',20*pi);
l.track(1,4).name = 'Rx4';                              % Set the MT1 name

l.track(1,5) = qd_track('circular',20*pi);
l.track(1,5).name = 'Rx5';                              % Set the MT1 name

l.track(1,6) = qd_track('circular',20*pi);
l.track(1,6).name = 'Rx6';                              % Set the MT1 name

l.track(1,7) = qd_track('circular',20*pi);
l.track(1,7).name = 'Rx7';                              % Set the MT1 name

l.track(1,8) = qd_track('circular',20*pi);
l.track(1,8).name = 'Rx8';                              % Set the MT1 name

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
% rx_loc = l.rx_position;
% plot scenario
if(show)
    l.visualize([],[],0); % Plot the layout
    view(-33, 60); % Enable 3D view
end

%% Set the scenario
l.set_scenario(par.scenario);
%% Generate channel coefficients
c = l.get_channels; % Generate channels
% cn = merge(c);
% cb = l.init_builder;                                    % Initialize channel builder object
% gen_ssf_parameters( cb );                               % Generate small-scale-fading parameters
% c = get_channels( cb );                                 % Get channel coefficients
% cn = merge(c2);
% c(8,1).coeff = c(8,1).coeff(:,:,:,1);
c_squeeze = zeros(par.B,c(1,1).no_snap,par.U);
for i = 1:par.U
    c_squeeze(:,:,i) = squeeze(c(i,1).coeff);
end
    
csi_mat = ones(par.B,par.U,size(c_squeeze,2));
for i = 1:size(csi_mat,3)
%      coll = permute(c(i,1).coeff,[2,1,3,4]);
%      coll = sum(coll,3);
     coll = permute(c_squeeze,[1,3,2]);
     csi_mat(:,:,i) = coll(:,:,i);
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
