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

This code should place the BS antenna array along the y-z axis as a plane facing East. The **par.array_v** parameter designates how many elements high the antenna array will be, and **par.array_h** is determined as the quotiet between the total number of antenna elements divided by **par.array_v**. For 128 antenna elements, the arrays I focussed on were 1x128, 2x64, 4x32, and 8x16.

The users are still being separated by a minimum angle of separation (described in other md file), are static, and are all 1.5m high. 

### All of the following simulations were performed at a carrier frequency of 3.5 GHz, bandwidth 10 MHz, 16QAM, and 8 users:

Images to be added in the morning.
