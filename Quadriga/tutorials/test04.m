%% Drifting Phases and Delays
%
% Drifting is the method used for obtaining time evolution within one segment. This tutorial
% demonstrates the effect of â€œdriftingâ€? on the channel coefficients. It shows how drifting can be
% enabled and disabled as well as how the resulting data can be analyzed.   
%
% Drifting is an essential feature of the channel model. Drifting enables a continuous time
% evolution of the path delays, the path phases, the departure- and arrival angles and the LSPs. It
% is thus the enabling feature for time continuous channel simulations. Although drifting was
% already available in the SCME branch of the WINNER channel model, it did not make it into the main
% branch. Thus, drifting is not available in the WIM1, WIM2 or WIM+ model. Here the functionality is
% implemented again. This script focuses on the delay and the phase component of the drifting
% functionality. 

%% Channel model setup and coefficient generation
% First, we parameterize the channel model. We start with the basic simulation parameters. For the
% desired output, we need two additional options: we want to evaluate absolute delays and we need to
% get all 20 sub-paths. Normally, the sub-paths are added already in the channel builder.   

clear all
close all

set(0,'defaultTextFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontName','Times')                    % Default Font Type
set(0,'defaultTextFontName','Times')                    % Default Font Type
set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type
set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Default Paper Size

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


%% Results and discussion
% The following plots represent the results of the test. The first plot shows the delay of the LOS
% tap (blue) and the delay of the first NLOS tap (red) vs. distance. The solid lines are from the
% channel with drifting, the dashed lines are from the channel without. The LOS delay is always
% increasing since the Rx is moving away from the Tx. However, the increase is not linear due to the
% 25 m height of the Tx. Without drifting, the delays are not updated and stay constant during the
% segment. The position of the first scatterer is in close distance to the Rx (only some m away).
% When moving, the Rx first approaches the scatterer (delay gets a bit smaller) and then the
% distance increases again.

set(0,'DefaultFigurePaperSize',[14.5 4.5])              % Change Paper Size
figure('Position',[ 100 , 100 , 760 , 400]);            % New figure

distance = c.rx_position(1,:);                          % 2D distance between Tx and Rx
plot( distance, c.delay(1,:)*1e9 , '-b' )               % Plot LOS delay with drifting
hold on
plot( distance, d.delay(1,:)*1e9 , '-.b' )              % Plot LOS delay without drifting
plot( distance, c.delay(2,:)*1e9 , '-r' )               % Plot 1st NLOS path with drifting
plot( distance, d.delay(2,:)*1e9 , '-.r' )              % Plot 1st NLOS path without drifting
hold off
xlabel('Distance from track start point')
ylabel('Delay [ns] ')
title('Path delays')
legend('LOS with drifting','LOS without drifting','NLOS with drifting','NLOS without drifting')

%%
% This plot shows the power of the first NLOS tap along the track. The fading is significantly
% higher in the beginning and becomes much less strong towards the end.

figure('Position',[ 100 , 100 , 760 , 400]);            % New figure
pow = abs(squeeze(sum( c.coeff(1,1,2,:,:) , 5 ))).^2;   % Calculate power of first NLOS path
plot( distance,10*log10(pow),'-r' )                     % Plot power of first NLOS path
xlabel('Distance from track start point')
ylabel('Tap power (dB)')
title('NLOS power with drifting')

%%
% Without drifting, the phases of the subpaths are approximated by assuming that the angles to the
% LBSs do not change. However, this only holds when the distance to the LBS is large. Here, the
% initial distance is small (ca. 5 m). When the initial angles are kept fixed along the track, the 
% error is significant. Here, the phase ramp is negative, indicating a movement direction towards
% the scatterer and thus a higher Doppler frequency. However, when the scatterer is passed, the Rx
% moves away from the scatterer and the Doppler frequency becomes lower. This is not reflected when
% drifting is turned off. 
%
% Note here, that with shorter delay spreads (as e.g. in satellite channels), the scatterers are
% placed closer to the Rxs initial position. This will amplify this effect. Hence, for correct time
% evolution results, drifting needs to be turned on.

figure('Position',[ 100 , 100 , 760 , 400]);            % New figure
pow = abs(squeeze(sum( d.coeff(1,1,2,:,:) , 5 ))).^2;   % Calculate power of first NLOS path
plot( distance,10*log10(pow),'-r' )                     % Plot power of first NLOS path
xlabel('Distance from track start point')
ylabel('Tap power (dB)')
title('NLOS power without drifting')

