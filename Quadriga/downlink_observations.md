Editted the following code to the downlink simulation Kaipeng provided me. Added the 'quadriga' case for the channel state matrix

```Matlab
    % generate iid Gaussian channel matrix and noise matrix
    N = sqrt(0.5)*(randn(par.U,par.T)+1i*randn(par.U,par.T));
    
    switch par.channel
        case 'rayleigh'
            H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));
        case 'los'
            [H_swm, H_pwm] = los(par);
            H = H_swm/norm(H_swm,'fro')*sqrt(par.B*par.U);
            %H = H_pwm/norm(H_pwm,'fro')*sqrt(par.B*par.U);
        case 'cellfree'
            [H] = cellfree(par);      
        case 'quadriga'
            H = channel_sim(par);
            norm_coef = zeros(1,par.U);                    
            for i = 1:par.U
                for j = 1:par.B
                    norm_coef(i)=norm_coef(i)+norm(H(j,i)); %sum the 2-norms of each column
                end
                norm_coef(i) = norm_coef(i)/par.B;         %average the 2-norm sum
                H(:,i) = H(:,i)/norm_coef(i);               %divide each entry of QuaDRiGa channel by avg 2-norm
            end
            H = H';
    end
 ```
channel_sim is the same generator that has been used in the simpleMIMOsim and decentralized simulator. The coefficients are normalized in the same way that was described in [this update.](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/observations_asof_July9.md) Since this simulator is on downlink, the channel coefficient matrix needs to be transposed.

Like with the uplink, the carrier frequency through which Quadriga generates the channel coefficients affects the effectiveness of the different downlink algorithms.

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc10e6.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc20e6.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc30e6.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc10e9.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc10e6.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc20e6.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc60e6.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc30e9.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/mmMagicLOS_fc10e6.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/mmMagicLOS_fc10e9.png" width="400" height="300">

## Also went on to test the difference between the encoded par.mod types.

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaLos_fc20e6_bpsk.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaLos_fc20e6_qpsk.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaLos_fc20e6_8psk.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaLos_fc20e6_16qam.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaLos_fc20e6_64qam.png" width="400" height="300">

## Also, there are some interesting effects from altering the number of decentralized divisions 

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc20e6_16qam_parC2.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc20e6_16qam_parC4.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc20e6_16qam_parC8.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc20e6_16qam_parC16.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc20e6_16qam_parC2.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc20e6_16qam_parC4.png" width="400" height="300">

<img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc20e6_16qam_parC8.png" width="400" height="300"> <img src="https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc20e6_16qam_parC16.png" width="400" height="300">

For the quadriga_MRT algorithm to start being more viable/close to the rayleigh approximation, the carrier frequency needs to approach the gigahertz regime. In the case of the mmMAGIC scenario, the quadriga_MRT algorithm actually surpasses the rayleigh_MRT in the gigahertz regime. The par.mod decisions yield expected results -> worse BER as the par.mod gets more complicated. The par.C parameter also yields results that are expected/found in previous simulator: the fewer the partitions, the better the results. The partitions were tested on a carrier frequency of 20 MHz, so if the carrier frequency was raised, there would be better results, but I thought this frequency would give a good demonstration/movement of the BER curves.
