Code that allows for dimensional configuration of antenna arrays:
```Matlab
% place BS antennas as a 2D array
BS_x_locs = zeros(par.array_v,par.array_h);
BS_y_locs = s.wavelength/2*(-(par.array_h-1)/2:1:(par.array_h-1)/2);
BS_z_locs = (s.wavelength/2*(-(par.array_v-1)/2:1:(par.array_v-1)/2))';

BS_y_locs = kron(ones(par.array_v,1),BS_y_locs);
BS_z_locs = 25*ones(par.array_v, par.array_h) + kron(ones(1,par.array_h),BS_z_locs);

BS_x_locs = reshape(BS_x_locs',par.B,1);
BS_y_locs = reshape(BS_y_locs',par.B,1);
BS_z_locs = reshape(BS_z_locs',par.B,1);
```
The results of this code can be seen below:

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/2x64_antenna.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/4x32_antenna.png" width="400" height="300">
<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/8x16_antenna.png" width="400" height="300">

This code should place the BS antenna array along the y-z axis as a plane facing East. The **par.array_v** parameter designates how many elements high the antenna array will be, and **par.array_h** is determined as the quotiet between the total number of antenna elements divided by **par.array_v**. For 128 antenna elements, the arrays I focussed on were 1x128, 2x64, 4x32, and 8x16.

The users are still being separated by a minimum angle of separation (described in other md file), are static, and are all 1.5m high. 

### All of the following simulations were performed at a carrier frequency of 3.5 GHz, bandwidth 10 MHz, 16QAM, and 8 users:
Images are in the order of 1x128, 2x64, 4x32, 8x16

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/1x128_8users.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/2x64_8users.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/4x32_8users.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/8x16_8users.png" width="400" height="300"> 

Taking a closer look at the effects of dimension across the algorithms:

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/fd_wf_comparison.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/mrt_comparison.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/wf_comparison.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/zf_comparison.png" width="400" height="300">

Then also with the affects of clusters:

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/pd_wf_comparison.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/pd_wf_comparison_parC4.png" width="400" height="300"> 

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/antenna%20geometries/pd_wf_comparison_parC8.png" width="400" height="300"> 
