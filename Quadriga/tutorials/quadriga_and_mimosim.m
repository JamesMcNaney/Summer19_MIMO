% generate par parameters to be used in simpleMIMOsim and channel_sim

par.simName = 'ERR_8x128_16QAM'; % simulation name (used for saving results)
par.runId = 0; % simulation ID (used to reproduce results)
par.MR = 128; % receive antennas 
par.MT = 8; % transmit antennas (set not larger than MR!) 
par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.trials = 250; % number of Monte-Carlo trials (transmissions)
par.SNRdB_list = -10:2:25; % list of SNR [dB] values to be simulated
% par.detector = {'ZF','bMMSE','uMMSE','ML'}; % define detector(s) to be simulated  
par.detector = {'ZF','bMMSE','uMMSE'};
% Have not been able to successfully run ML

% To be used in channel_sim.m -> where QuaDRiGa generates coefficients
par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
par.fc = 28e6; % carrier frequency [Hz]
par.BW = 10e6; % bandwidth [Hz]
par.N = 1024; % number of carriers
par.B = par.MR; % number of antennas in the BS (we use a single BS)
par.U = par.MT; % number of single-antenna UEs
par.iid = 1;    % simulates Gaussian iid

% Run the simulation

% subplot(1,2,1)
simpleMIMOsim(par); hold on;
% title('Gaussian iid channel');

% subplot(1,2,2)
par.iid = 0;    % simulates the par.scenario
simpleMIMOsim(par);
hold off;
title(['Gaussian iid and' string(par.scenario)], 'Interpreter', 'none');
% title('QuaDRiGa BERLIN\_UMa\_NLOS Channel');

    