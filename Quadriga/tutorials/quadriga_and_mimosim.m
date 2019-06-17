% generate par parameters to be used in simpleMIMOsim and channel_sim

    par.simName = 'ERR_MTxMR_16QAM'; % simulation name (used for saving results)
    par.runId = 0; % simulation ID (used to reproduce results)
    par.MR = 64; % receive antennas 
    par.MT = 8; % transmit antennas (set not larger than MR!) 
    par.mod = 'QPSK'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 500; % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = -10:2:25; % list of SNR [dB] values to be simulated
%     par.detector = {'ZF','bMMSE','uMMSE','ML'}; % define detector(s) to be simulated  
    par.detector = {'bMMSE','ZF'};

par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
par.fc = 60e9; % carrier frequency [Hz]
par.BW = 14e6; % bandwidth [Hz]
par.N = 2048; % number of carriers
par.B = 64; % number of antennas in the BS (we use a single BS)
par.U = 8; % number of single-antenna UEs
par.iid = 1;
% Run the simulation

subplot(1,2,1)
simpleMIMOsim(par);
title('Gaussian iid channel');

subplot(1,2,2)
par.iid = 0;
simpleMIMOsim(par);
title('QuaDRiGa Freespace Channel');

    