s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 2.53e9;                            % 2.53 GHz carrier frequency
s.sample_density = 4;                                   % 4 samples per half-wavelength
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.show_progress_bars = 0;                               % Disable progress bars

l = qd_layout(s);
l.no_rx = 10;
l.tx_array = qd_arrayant('dipole');                     % Dipole antennas at all Rx and Tx
l.rx_array = l.tx_array;
l.tx_position(3) = 25;  

l.track(1,1) = qd_track('circular',10*pi,0);
l.track(1,2) = qd_track('circular',10*pi,0);
t.initial_position = [5;0;0];
t.segment_index = [1,40,90];
t.scenario = { 'C2l' , 'C2n' , 'C2l' };
t.interpolate_positions( 100 );
t.compute_directions;