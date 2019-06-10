function newf = proof_of_concept()

gridx = 7;
gridy = 7;
gridz = 1;
rx_space = 2;


set(0,'defaultTextFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontName','Times')                    % Default Font Type
set(0,'defaultTextFontName','Times')                    % Default Font Type
set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type

s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 2.53e9;                            % 2.53 GHz carrier frequency
s.sample_density = 4;                                   % 4 samples per half-wavelength
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.show_progress_bars = 0;                               % Disable progress bars

%%
% a = qd_arrayant('lhcp-rhcp-dipole');    % Create circular polarized antenna
% 
% a2 = qd_arrayant('custom',90,90,0);     % Create linear polarized patch antenna
% a2.copy_element(1,2);                   % Copy the antenna element
% a2.rotate_pattern(90,'x',2);            % Rotate second element by 90 degree
% 
% a.append_array( a2 );                   % Append the second antenna to the first

a = qd_arrayant('multi', 8, 0.5, 12 );
a.no_elements = 8;
a2 = qd_arrayant('omni', 1);

% a.visualize;
% a2.visualize;

%%
% Second, we create a more complex network layout featuring an elevated transmitter (25 m) and two
% receivers at 1.5 m height. The first Rx moves along a circular track around the receiver. The
% second receiver moves away from the Tx. Both start at the same point. Note here, that each track
% is split into three segments. The first Rx goes from an LOS area to a shaded area and back. The
% second track also start in the LOS area. Here, the scenario changes to another LOS segment and
% then to an NLOS segment. The LOS-LOS change will create new small-scale fading parameters, but the
% large scale parameters (LSPs) will be highly correlated between those two segments.

dummy = gridx*gridy;
l = qd_layout(s);                                       % Create new QuaDRiGa layout
l.no_rx = gridx*gridy*gridz;                    %james: Initialize set number of receivers (made 3d array of receivers for calculation)        
l.tx_array = a;                                 %james: these are the two lines that make all channel coefficients 4d complex       
l.rx_array = a2;                                 %james: but break the pdp calculations
% l.tx_array = qd_arrayant('dipole');           % Dipole antennas at all Rx and Tx
% l.rx_array = l.tx_array;
l.tx_position(3) = 25;                          % Elevate Tx to 25 m
l.track = qd_track('linear',0,0);


%UMal = 'BERLIN_UMa_LOS';                                % LOS scenario name
%UMan = 'BERLIN_UMa_NLOS';                               % NLOS scenario name
UMal = 'LOSonly';

%%
iter = 1;
gridx2 = gridx/2;
gridy2 = gridy/2;
gridz2 = gridz/2;
for rx_x = -gridx2:gridx2
    for rx_y = -gridy2:gridy2
        for rx_z = 1:gridz
            l.rx_position(1,iter)  = rx_x*rx_space;
            l.rx_position(2,iter)  = rx_y*rx_space;
            l.rx_position(3,iter)  = rx_z*rx_space;
            iter = iter + 1;
        end    
    end
end
l.set_scenario(UMal);                      % Use UMal scenario
%%
set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Adjust paper size for plot
% l.visualize;                                            % Plot the layout

% compute_directions( l.track );                          % Align antenna direction with track

%%
% Now we create the channel coefficients. The fixing the random seed guarantees repeatable results
% (i.e. the taps will be at the same positions for both runs).

p = l.init_builder;                                     % Create channel builders


%%

p.gen_ssf_parameters;                                   % Generate small-scale fading

%%

%s.use_spherical_waves = 1;                              % Enable drifting (=spherical waves)
dist = get_distances( p );                               

%d = l.get_channels;                                     %what are the differences between this and the next line?
c = get_channels( p );                                   %because the channel coefficients differ when generated from layout vs qd_builder
los_c = get_los_channels( p );              %since all receivers are LOS, this merges all coefficients together

cn = merge( c );

freq_response = cn(1,1).fr(100e6,64);
newf = reshape(freq_response,64,8);
end