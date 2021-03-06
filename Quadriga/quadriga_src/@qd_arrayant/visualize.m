function h_figures = visualize( h_qd_arrayant, i_element )
%VISUALIZE Create a plot showing the element configurations
%
% Calling object:
%   Single object
%
% Input:
%   i_element
%   The element indices for which the plot os created. If no element index are given, a plot is
%   created for each element in the array. 
%
% Output:
%   h_figures
%   The figure handles for further processing of the images.
%
%
% QuaDRiGa Copyright (C) 2011-2017 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if numel( h_qd_arrayant ) > 1 
   error('QuaDRiGa:qd_arrayant:visualize','visualize not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

if nargin == 1
    i_element = 1:h_qd_arrayant.no_elements;
end

if ~(any(size(i_element) == 1) && isnumeric(i_element) ...
        && isreal(i_element) && all(mod(i_element, 1) == 0) && all(i_element > 0))
    error('??? "element" must be integer and > 0')
elseif any(i_element > h_qd_arrayant.no_elements)
    error('??? "element" exceed "no_elements"')
end

az = double(h_qd_arrayant.azimuth_grid);
if numel(az) == 361
    indt = 1:5:361;
else
    indt = 1:numel(az);
end
az = az(indt);

elev = double(h_qd_arrayant.elevation_grid');
if numel(elev) == 181
    indp = 1:5:181;
else
    indp = 1:numel(elev);
end
elev = elev(indp);

[az_grid, elev_grid] = meshgrid(az, elev);
[Xi, Yi, Zi] = sph2cart(az_grid, elev_grid, 1);

min_value = -20;
scaling = 1;

h_figures = zeros(1, numel(i_element));

for n = 1:numel(i_element)
    
    % Read the element patterns
    Fa = h_qd_arrayant.Fa(indp, indt, i_element(n));
    Fb = h_qd_arrayant.Fb(indp, indt, i_element(n));
    
    % calculate radiation power pattern
    P = abs(Fa).^2 + abs(Fb).^2;
    
    % Normalize by max value
    P_max = max( P(:) );
    P = P ./ P_max;
    
    % Calculate directivity
    tmp         = cos(elev_grid(:));
    directivity_lin    = sum(tmp) ./ sum(P(:).*tmp);
    directivity_dbi    = 10*log10(directivity_lin);
    
    % Normalize patterns by directivity
    tmp = sqrt(directivity_lin./P_max);
    Fa = Fa .* tmp;
    Fb = Fb .* tmp;
    
    h_figures(n) = figure('Position', [50, 400, 1200, 500],...
        'Name', [h_qd_arrayant.name ' element ', num2str(i_element(n))]);
    
    axes('position',[0 0 0.92 0.9]);%,'Visible','Off');
    axis off
    
    title(['Array Antenna Element ', num2str(i_element(n))] );
    
    text(0.02,1,'D^{[\theta]}(\theta, \phi)');
    text(0.9 ,1,'D^{[\phi]}(\theta, \phi)');
    no_plots = 2;
    
    
    for m = 1:no_plots
        if no_plots == 3
            axes('position',[-0.25+0.3*m 0.12 0.25 0.7]);
            Po = 10*log10(abs(Fa).^2);
        else
            axes('position',[-0.4+0.45*m 0.12 0.38 0.82]);
            Po = 10*log10(abs(Fa).^2);
        end
        
        switch m
            case 1
                Po = 10*log10( abs(Fa).^2 );
            case 2
                Po = 10*log10( abs(Fb).^2 );
            case 3
                Po = 10*log10( abs(Fc).^2 );
        end
        Po = double( Po );
        
        P = Po;
        P(P < min_value) = min_value;
        P = (P - min_value) ./ (directivity_dbi - min_value) .* scaling;
        
        X = P .* Xi;
        Y = P .* Yi;
        Z = P .* Zi;
        
        surf(X, Y, Z, Po)
        
        axis equal
        axis(scaling.*[-1 1 -1 1 -1 1]);
        caxis([min_value, directivity_dbi]);
        set(gca, 'xtick', (-1:1).*scaling/2);
        set(gca, 'ytick', (-1:1).*scaling/2);
        set(gca, 'ztick', (-1:1).*scaling/2);
        xlabel('x')
        ylabel('y')
        zlabel('z')
       
        view(45, 33)
    end
    
    axes('position', [0.08 0.08 0.08 0.08],'Visible','Off');
    caxis([min_value, directivity_dbi]);
    han = colorbar('EastOutside', 'XTick', min_value:3:floor(directivity_dbi));
    set(han, 'position', [0.92 0.06 0.02 0.90])
    zlab = get(han, 'ylabel');
    set(zlab, 'String', 'Partial Directivity in dBi');
    
end


end
