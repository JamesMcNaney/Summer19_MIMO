# This document is an attempt to record my understanding of the **simpleMIMOsim** function Kaipeng provided. It is also my first attempt at becoming more versed in utilizing markdown.

Code being used: 
* [channel_sim.m](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/channel_sim.m)
* [simpleMIMOsim.m](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/simpleMIMOsim.m)
* [quadriga_and_mimosim.m](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/quadriga_and_mimosim.m)

As of June 14, 2019, I am utilizing this function in its default state while trying to integrate channel coefficients derived from the QuaDRiGa program. Therefore, the default simulation parameters are used and replicated in the QuaDRiGa setup. These parameters can be manipulated by the user by altering the corresponding values found within specific properties of the **par** struct. Of note: MR corresponds to the number of elements on the antenna array, MT corresponds to mobile terminals (users), mod corresponds to the modulation type -> alphabet used to transmit symbols, trails to the number of trials, SNRdb_list to the array of SNR values to be tested against, and the detector method to be used. This program is mainly meant to explore the Bit Error Rate (BER) as a function of the SNR and the output of this function (currently) is a plot of this relationship.
Further non-user-dependent initialization preceeds the *start simulation* stage of the program. A random bitstream is generated such that each trial will have 1 symbol per MT (number of bits being determined by the length needed for the modulation scheme). 
For each trial:
* The the index value (**idx**) that the random bitstream corresponds to is definitively expressed such that the symbol **s** which the bitstream should be recognized can be accessed from the **par**'s symbol property. Each entry of **s** (one for each MT) corresponds to a symbol sent from the MT to the BS (uplink).
* The noise and channel matrix are next with the noise matrix being a gaussian vector. *This is where QuaDRiGa's output can be implemented.* In the default state, **H** is generated as an iid Gaussian matrix (with a factor of 1/sqrt(2) in front... Normalization factor?). **H** has one entry per each MT-antenna element pair, which leads me to believe that the coefficients are just the complex *gain* coefficients between the two... Therefore, the values within the CSI matrix to be returned from QuaDRiGa should probably be the time-dependent gain coefficients that are generated initially and do not require further post-processing steps (like taking into the frequency domain).
* The ideal signal **x** that gets transmitted is calculated by multiplying H\*s. When initially transferring the coefficients from QuaDRiGa, all of the values were on the micro scale, which would kill the ideal transmit signal/have the noise cancel everything out. That is why I include a normalization of the CSI matrix (**I do not know if this is valid**), somewhat like how the iid Gaussian included a 1/sqrt(2) normalization factor. 
* Given the ideal signal **x**, it is time to iterate over the list of SNR values the user defines in the initial parameterization step. For each SNR value:
  * calculate the noise variance **N0**. 
  * calculate what the signal on uplink will look like to the BS: y = x+sqrt(N0)\*n
  * then identify which detector is being used (zero-forcing etc) and run the respective function for the detector. Store idxhat and bithat (to be used later)
  * determine whether idxhat is the same as the predetermined **idx** from earlier in the simulation loop
  * calculate respective error rates and store them in their respective properties within **res** struct
* Plot the BER vs SNR curve. 
### Following are a few of these curves using QuaDRiGa's (normalized) channel coefficients. Changes of the figures come mostly from reducing number of trials and SNR values tested across. There are other figures available in the tutorials/figures_images folder. These plots do not approach zero as quickly as the iid Gaussian estimate.

![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/quadriga_in_simulation.png)

![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/quadriga_in_simulation_2.png)

![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/quadriga_in_simulation_3.png)
---
![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/comparing_iidQuadriga_64x8.png)
Testing/comparing QPSK
![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/comparing_iidQuadriga_64x8_16QAM.png)
Testing/comparing 16QAM

Now, when testing different scenarios (I assume 'Freespace' means completely line of sight), I was running into dimension problems with the channel state matrix. This is because QuaDRiGa returns a coefficient for each generated path. There were a total of 25 paths generated with the Berlin NLOS scenario, and these different paths correspond to the paths the signal takes off of multiple scatterers. Since NLOS is one of the strengths of massive MUMIMO, I did not want to simply *choose* a single path's coefficients, so I summed them together. (**I do not know if this is nonsensical or not since I am not sure what the algorithms are supposed to do to multipath signals at the BS on UL**). The following was the result:
![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/Berlin_NLOS_summed.png)

Following are more results from the summed coefficients along the multipath dimension. *If* this process is correct, it's promising
![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/summed_freespace.png)
![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/summed_mmMagic_LOS.png)
![alt text](https://github.com/JamesMcNaney/Summer19_MIMO/blob/master/Quadriga/tutorials/figures_images/summed_mmMagic_NLOS.png)




## Plots and figures for some simulations can be found in the tutorials/figures_images folder.
