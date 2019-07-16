par.B = 128;
par.fc = 3.5e9;
par.array_v = 1;        % number of vertical antenna elements, assuming total
                    % number of antenna elements is a power of 2, must also
                    % be a factor of 2 < total number of antenna elements

par.array_h = par.B/par.array_v;

s = qd_simulation_parameters;
s.sample_density = 2; % 2 samples per half-wavelength (spacial sampling)
s.center_frequency = par.fc;
s.use_absolute_delays = false; % Include delay of the LOS path
s.use_spherical_waves = true; % No need to use spherical waves if rx and tx are far enough
s.show_progress_bars = false;

% original antenna placement code
% BS_x_locs = zeros(par.B,1);
% BS_y_locs = s.wavelength/2*(-(par.B-1)/2:1:(par.B-1)/2)';
% BS_z_locs = 25*ones(par.B,1); % note that the height of the transmitter is 25m

% new antenna placement code
BS_x_locs = zeros(par.array_v,par.array_h);
BS_y_locs = s.wavelength/2*(-(par.array_h-1)/2:1:(par.array_h-1)/2);
BS_z_locs = (s.wavelength/2*(-(par.array_v-1)/2:1:(par.array_v-1)/2))';

BS_y_locs = kron(ones(par.array_v,1),BS_y_locs);
BS_z_locs = kron(ones(1,par.array_h),BS_z_locs);

BS_x_locs = reshape(BS_x_locs',128,1);
BS_y_locs = reshape(BS_y_locs',128,1);
BS_z_locs = reshape(BS_z_locs',128,1);

scatter3(BS_x_locs, BS_y_locs, BS_z_locs);