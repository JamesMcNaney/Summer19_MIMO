%% Effects of the Antenna-Orientation

%% Model and Antenna Setup
% Here, we parametrize the simulation. We place the receiver 10 m away from the transmitter and
% chose the scenario "LOSonly". Thus, no NLOS components are present. The receiver is set up as a
% multi-element array antenna using both, circular and linear polarization.

clear all
close all

s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 30e9;                              % 30 GHz carrier frequency (mmWave regime)
s.sample_density = 2.5;                                 % 2.5 samples per half-wavelength
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.show_progress_bars = 0;                               % Disable progress bars

a = qd_arrayant('3gpp-3d',1,64,30e9,1,pi/4,.5,1,1);    % Create 1x64 3gpp-mmw antenna array
%%
t = qd_track('linear',20,0);                        % 20 m track, direction SE
t.initial_position = [0;0;1.5];                        % Start position
t.interpolate_positions( 128/20 );                      % Interpolate
% t.segment_index       = [1,40,90];                      % Assign segments
t.scenario            = '3GPP_38.901_UMi_LOS';
t.interpolate_positions( s.samples_per_meter );         % Apply sample density

l = qd_layout( s );                                     % New QuaDRiGa layout
l.tx_array = a;                                         % Set 3gpp-mmw antenna
l.tx_array.rotate_pattern(90,'z');

l.rx_array = qd_arrayant('omni');                       % Set Dipole antenna
l.tx_position = [10;-5;10];
l.track = t;                                            % Assign track

set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Adjust paper size for plot
l.visualize;                                            % Plot the layout

%%
% l = qd_layout;
% l.simpar.show_progress_bars = 0;        % Disable progress bar indicator
% 
% t = qd_track('linear',20,0);      % Create new track (pi turns the rx by 180 degree)
% t.initial_position = [0;0;0];               % Set the receiver position
% t.interpolate_positions( 128/20 );
% t.scenario = '3GPP_38.901_UMi_LOS';
% 
% l.tx_position = [10;-5;10];
% 
% l.track = t;
% 
% l.tx_array = a;             
% l.rx_array = qd_arrayant('omni');
% l.tx_array.rotate_pattern(90,'z');
% % l.track.no_snapshots = 20;
% t.set_speed(20);
% l.visualize; 
% test = l.get_channels;
%%
cn = l.get_channels;                                    % Generate channels

t.set_speed( 30 );                                       % Set constant speed
dist = t.interpolate_movement( s.wavelength/(2*20) );   % Get snapshot positions
ci = cn.interpolate( dist , 'spline' );                 % Interpolate channels
