Code that allows for dimensional configuration of antenna arrays:
```Matlab
% place BS antennas as a 2D array
BS_x_locs = zeros(par.array_v,par.array_h);
BS_y_locs = s.wavelength/2*(-(par.array_h-1)/2:1:(par.array_h-1)/2);
BS_z_locs = (s.wavelength/2*(-(par.array_v-1)/2:1:(par.array_v-1)/2))';

BS_y_locs = kron(ones(par.array_v,1),BS_y_locs);
BS_z_locs = 25*ones(par.array_v, par.array_h) + kron(ones(1,par.array_h),BS_z_locs);

BS_x_locs = reshape(BS_x_locs',128,1);
BS_y_locs = reshape(BS_y_locs',128,1);
BS_z_locs = reshape(BS_z_locs',128,1);
```
Images to be added in the morning.
