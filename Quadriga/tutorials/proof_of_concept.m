function newf = proof_of_concept()

% gridx = 7;
% gridy = 7;
% gridz = 1;
% rx_space = 2;
% 
% 
% set(0,'defaultTextFontSize', 18)                        % Default Font Size
% set(0,'defaultAxesFontSize', 18)                        % Default Font Size
% set(0,'defaultAxesFontName','Times')                    % Default Font Type
% set(0,'defaultTextFontName','Times')                    % Default Font Type
% set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
% set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type
% 
% s = qd_simulation_parameters;                           % New simulation parameters
% s.center_frequency = 2.53e9;                            % 2.53 GHz carrier frequency
% s.sample_density = 4;                                   % 4 samples per half-wavelength
% s.use_absolute_delays = 1;                              % Include delay of the LOS path
% s.show_progress_bars = 0;                               % Disable progress bars
% 
% %%
% % a = qd_arrayant('lhcp-rhcp-dipole');    % Create circular polarized antenna
% % 
% % a2 = qd_arrayant('custom',90,90,0);     % Create linear polarized patch antenna
% % a2.copy_element(1,2);                   % Copy the antenna element
% % a2.rotate_pattern(90,'x',2);            % Rotate second element by 90 degree
% % 
% % a.append_array( a2 );                   % Append the second antenna to the first
% 
% a = qd_arrayant('dipole', 8, 0.5, 12 );
% a.no_elements = 8;
% a2 = qd_arrayant('omni', 1);
% 
% % a.visualize;
% % a2.visualize;
% 
% %%
% % Second, we create a more complex network layout featuring an elevated transmitter (25 m) and two
% % receivers at 1.5 m height. The first Rx moves along a circular track around the receiver. The
% % second receiver moves away from the Tx. Both start at the same point. Note here, that each track
% % is split into three segments. The first Rx goes from an LOS area to a shaded area and back. The
% % second track also start in the LOS area. Here, the scenario changes to another LOS segment and
% % then to an NLOS segment. The LOS-LOS change will create new small-scale fading parameters, but the
% % large scale parameters (LSPs) will be highly correlated between those two segments.
% 
% dummy = gridx*gridy;
% l = qd_layout(s);                                       % Create new QuaDRiGa layout
% l.no_rx = gridx*gridy*gridz;                    %james: Initialize set number of receivers (made 3d array of receivers for calculation)        
% l.tx_array = a;                                 %james: these are the two lines that make all channel coefficients 4d complex       
% l.rx_array = a2;                                 %james: but break the pdp calculations
% % l.tx_array = qd_arrayant('dipole');           % Dipole antennas at all Rx and Tx
% % l.rx_array = l.tx_array;
% l.tx_position(3) = 25;                          % Elevate Tx to 25 m
% l.track = qd_track('linear',0,0);
% 
% 
% %UMal = 'BERLIN_UMa_LOS';                                % LOS scenario name
% %UMan = 'BERLIN_UMa_NLOS';                               % NLOS scenario name
% UMal = 'LOSonly';
% 
% %%
% iter = 1;
% gridx2 = gridx/2;
% gridy2 = gridy/2;
% gridz2 = gridz/2;
% for rx_x = -gridx2:gridx2
%     for rx_y = -gridy2:gridy2
%         for rx_z = 1:gridz
%             l.rx_position(1,iter)  = rx_x*rx_space;
%             l.rx_position(2,iter)  = rx_y*rx_space;
%             l.rx_position(3,iter)  = rx_z*rx_space;
%             iter = iter + 1;
%         end    
%     end
% end
% l.set_scenario(UMal);                      % Use UMal scenario
% %%
% set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Adjust paper size for plot
% % l.visualize;                                            % Plot the layout
% 
% % compute_directions( l.track );                          % Align antenna direction with track
% 
% %%
% % Now we create the channel coefficients. The fixing the random seed guarantees repeatable results
% % (i.e. the taps will be at the same positions for both runs).
% 
% p = l.init_builder;                                     % Create channel builders
% 
% 
% %%
% 
% p.gen_ssf_parameters;                                   % Generate small-scale fading
% 
% %%
% 
% %s.use_spherical_waves = 1;                              % Enable drifting (=spherical waves)
% dist = get_distances( p );                               
% 
% %d = l.get_channels;                                     %what are the differences between this and the next line?
% c = get_channels( p );                                   %because the channel coefficients differ when generated from layout vs qd_builder
% los_c = get_los_channels( p );              %since all receivers are LOS, this merges all coefficients together
% 
% los_c.coeff = (reshape(los_c.coeff,8,64))';
% 
% cn = merge( c );
% 
% freq_response = los_c.coeff;
% newf = freq_response;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 2.53e9;                            % 2.53 GHz carrier frequency
s.sample_density = 4;                                   % 4 samples per half-wavelength
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.show_progress_bars = 0;                               % Disable progress bars

l = qd_layout( s );                                     % Create a new Layout
l.tx_array = qd_arrayant('3gpp-3d',4,2);                     % create V-polarized dipole
% l.tx_array.no_elements = 8;
l.tx_array.set_grid( (-180:10:180)*pi/180 , (-90:10:90)*pi/180 );
l.tx_array.Fa = l.tx_array.Fa ./ max(l.tx_array.Fa(:));

% l.tx_array.copy_element(1,2:8);                         % Duplicate the elements
% l.tx_array.rotate_pattern(45,'y',2);                    % 45 degree polarization
% l.tx_array.rotate_pattern(90,'y',3);                    % 90 degree polarization
% l.rx_array = l.tx_array;                                % Use the same array for the Rx
l.rx_array = qd_arrayant('dipole');

set(0,'DefaultFigurePaperSize',[14.5 5.3])              % Adjust paper size for plot
% l.tx_array.visualize(1);pause(1);                       % Plot the first antenna element
% l.tx_array.visualize(2);pause(1);                       % Plot the second antenna element
% l.tx_array.visualize(3);pause(1);                       % Plot the third antenna element

%% Defining a track
% The third step defines the track. Here, we use a circle with 40 m diameter starting in the east,
% traveling north. We also choose a LOS scenario since we want to study the LOS polarization. The
% transmitter is located 12 m north of the center of the circle at an elevation of 6 m.

l.track = qd_track('circular',0,0);                 % Circular track, radius 20 m
% interpolate_positions( l.track, s.samples_per_meter );  % Interpolate positions
l.tx_position = [ 0 ; 12 ; 6 ];                         % Tx position
l.rx_position = [ 20 ; 0 ; 0 ];                         % Start position for the Rx track
l.set_scenario('BERLIN_UMa_LOS');

set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Adjust paper size for plot
l.visualize;                                            % Plot the layout

%%
% a = qd_arrayant('dipole', 8, 0.5, 12 );
% a.no_elements = 8;
% a2 = qd_arrayant('omni', 1);
% 
% a.visualize;
% a2.visualize;
% 
% l.tx_array = a;
% l.rx_array = a2;
%%
% Now, we generate the LSPs. We set the shadow fading and K-factor to 1 and disable the path loss
% model. 

cb = l.init_builder;                                    % Create new builder object
cb.scenpar.SF_sigma = 0;                                % 0 dB shadow fading
cb.scenpar.KF_mu = 0;                                   % 0 dB K-Factor
cb.scenpar.KF_sigma = 0;                                % No KF variation
cb.plpar = [];                                          % Disable path loss model
cb.gen_ssf_parameters;                                  % Generate large- and small-scale fading

%%
% Now, we generate the channel coefficients. The first run uses the drifting module, the second run
% disables it. Note that drifting needs significantly more computing resources. In some scenarios it
% might thus be useful to disable the feature to get quicker simulation results.   

s.use_spherical_waves = 1;                              % Enable drifting (=spherical waves)
c = cb.get_channels;                                    % Generate channel coefficients
c.individual_delays = 0;                                % Remove per-antenna delays

s.use_spherical_waves = 0;                              % Disable drifting
d = cb.get_channels;                                    % Generate channel coefficients

newf = zeros(64,8);
for i = 1:64
    newf(i,:) = c.coeff(:,:,i,1);
end

end