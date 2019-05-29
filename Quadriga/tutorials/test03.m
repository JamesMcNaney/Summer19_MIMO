%% Channel model setup and coefficient generation
% First, we set up the channel model.

close all
clear all

gridx = 10;
gridy = 10;
rx_space = 1;


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
% Second, we create a more complex network layout featuring an elevated transmitter (25 m) and two
% receivers at 1.5 m height. The first Rx moves along a circular track around the receiver. The
% second receiver moves away from the Tx. Both start at the same point. Note here, that each track
% is split into three segments. The first Rx goes from an LOS area to a shaded area and back. The
% second track also start in the LOS area. Here, the scenario changes to another LOS segment and
% then to an NLOS segment. The LOS-LOS change will create new small-scale fading parameters, but the
% large scale parameters (LSPs) will be highly correlated between those two segments.

dummy = gridx*gridy;
l = qd_layout(s);                                       % Create new QuaDRiGa layout
l.no_rx = gridx*gridy;                                            % Two receivers
l.tx_array = qd_arrayant('dipole');                     % Dipole antennas at all Rx and Tx
l.rx_array = l.tx_array;
l.tx_position(3) = 25;                                  % Elevate Tx to 25 m

UMal = 'BERLIN_UMa_LOS';                                % LOS scenario name
UMan = 'BERLIN_UMa_NLOS';                               % NLOS scenario name

iter = 1;
for rx_x = 1:gridx
    for ry_y = 1:gridy
%         l.track(1,iter) = qd_track('/iter' , 0, 0);                % linear track with 0m length... static
        l.rx_position(1,iter)  = rx_x;
        l.rx_position(2,iter)  = ry_y;
        l.rx_position(3,iter)  = 1.5;
%         l.track(1,iter).scenario          = UMal;                  % Scenario
        iter = iter + 1;
    end
end
l.set_scenario(UMal);                      % Use UMal scenario

set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Adjust paper size for plot
l.visualize;                                            % Plot the layout

% compute_directions( l.track );                          % Align antenna direction with track

%%
% Now we create the channel coefficients. The fixing the random seed guarantees repeatable results
% (i.e. the taps will be at the same positions for both runs).

% p = l.init_builder;                                     % Create channel builders
p = l.init_builder;

%%

p.gen_ssf_parameters;                                   % Generate small-scale fading

%%

s.use_spherical_waves = 1;                              % Enable drifting (=spherical waves)
d = l.get_channels;                                     %what are the differences between this and the next line?
c = get_channels( p ); 
e = p.get_channels;                                     %same as line above
cn = merge( c ); 

%%
% Next, we plot the power-delay profiles for both tracks. We calculate the frequency response of the
% channel and transform it back to time domain by an IFFT. Then, we create a 2D image of the
% received power at each position of the track. We start with the circular track.

h = cn(1,1).fr( 100e6,512 );                            % Freq.-domain channel
h = squeeze(h);                                         % Remove singleton dimensions
pdp = 10*log10(abs(ifft(h,[],1).').^2);                 % Power-delay profile

figure('Position',[ 100 , 100 , 760 , 400]);            % New figure
imagesc(pdp(:,1:256));

caxis([ max(max(pdp))-50 max(max(pdp))-5 ]); colorbar;  % Figure decorations
cm = colormap('hot'); colormap(cm(end:-1:1,:));
set(gca,'XTick',1:32:255); set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
set(gca,'YTick',1:cn(1).no_snap/8:cn(1).no_snap);
set(gca,'YTickLabel', (0:cn(1).no_snap/8:cn(1).no_snap)/cn(1).no_snap * 360 );
xlabel('Delay [\mus]'); 
title('PDP1');

%%
% The X-axis shows the delay in microseconds and the Y-axis shows the position on the circle. For
% easier navigation, the position is given in degrees. 0 deg means east (starting point), 90 deg
% means north, 180 deg west and 270 deg south. The LOS delay stays constant since the distance to
% the Tx is also constant. However, the power of the LOS changes according to the scenario. Also
% note, that the NLOS segment has more paths due to the longer delay spread.
%
% Next, we create the same plot for the linear track. Note the slight increase in the LOS delay and
% the high similarity of the first two LOS segments due to the correlated LSPs. Segment change is at
% around 6 m.

h = cn(1,2).fr( 100e6,512 );                            % Freq.-domain channel
h = squeeze(h);                                         % Remove singleton dimensions
pdp = 10*log10(abs(ifft(h,[],1).').^2);                 % Power-delay profile

figure('Position',[ 100 , 100 , 760 , 400]);            % New figure
imagesc(pdp(:,1:256));

caxis([ max(max(pdp))-50 max(max(pdp))-5 ]); colorbar;  % Figure decorations
cm = colormap('hot'); colormap(cm(end:-1:1,:));
set(gca,'XTick',1:32:255); set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
set(gca,'YTick',1:cn(2).no_snap/8:cn(2).no_snap);
set(gca,'YTickLabel', (0:cn(2).no_snap/8:cn(2).no_snap)/cn(2).no_snap * 20 );
xlabel('Delay [\mus]'); 
title('PDP2');