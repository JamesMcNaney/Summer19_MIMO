 par.runId = 0; % simulation ID (used to reproduce results)
    par.C = 4;
    par.MR = 128; % total receive antennas
    
    par.MT = 16; % transmit antennas (set not larger than MR!)
    
    par.beta = par.MT/par.MR;
    % number of clusters (uniformly distributed antennas over C)
    % (set C>1 as C=1 gives an error for the decentralized detectors
    % MR/C has to be integer!

    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 250; % number of Monte-Carlo trials (transmissions)    
    
%% added parameters for QuaDRiGa
    par.scenario = 'mmMAGIC_UMi_LOS'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
    par.fc = 3.5e9; % carrier frequency [Hz]
    par.BW = 10e6; % bandwidth [Hz]
    par.N = 1024; % number of carriers
    par.B = par.MR; % number of antennas in the BS (we use a single BS)
    par.U = par.MT; % number of single-antenna UEs
    par.normalize = 1;
%% sim parameters (please read!)
    par.SNRdB_list = [-15:5:20]; % list of SNR [dB] values to be simulated
    par.detector = {...         
         'uMMSE',...
         'uMMSE_decent',...
         'CG',...
         'DCG',...
         'MF',...
        }; % define detector(s) to be simulated
    
    
    % CG
    par.iteration = 4; % number of iteration for CG
    par.normalize = 1;
    par.simName = ['results/' num2str(par.MR) 'x' num2str(par.MT)...
        '_' num2str(par.C) 'clusters_' num2str(par.trials) 'trials'];
    
    
%%
par.iid = 0;
par.array_v = 1;
par.array_h = par.B/par.array_v;
decentralized_MIMOsim(par);
% hold on
% par.iid = 0;
% decentralized_MIMOsim(par);
% hold off
% title(['Par.C = ' num2str(par.C) par.scenario ], 'Interpreter','none');