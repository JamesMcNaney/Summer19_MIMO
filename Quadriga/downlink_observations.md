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

Like with the uplink, the carrier frequency through which Quadriga generates the channel coefficients affects the effectiveness of the different downlink algorithms.


![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc10e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc20e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc30e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/berlinUmaNlos_fc10e9.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc10e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc20e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc60e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/freespace_fc30e9.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/mmMagicLOS_fc10e6.png)

![alt_text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/downlink_figs/mmMagicLOS_fc10e9.png)
