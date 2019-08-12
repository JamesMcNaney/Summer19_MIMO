clear all 

par.runId = 0; % simulation ID (used to reproduce results)
    par.C = 4;
    par.MR = 128; % total receive antennas
    
    par.MT = 16; % transmit antennas (set not larger than MR!)
    
    par.beta = par.MT/par.MR;
    % number of clusters (uniformly distributed antennas over C)
    % (set C>1 as C=1 gives an error for the decentralized detectors
    % MR/C has to be integer!

    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 2000; % number of Monte-Carlo trials (transmissions)    
    
%% added parameters for QuaDRiGa
    par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
    par.fc = 3.5e9; % carrier frequency [Hz]
    par.BW = 10e6; % bandwidth [Hz]
    par.N = 1024; % number of carriers
    par.B = par.MR; % number of antennas in the BS (we use a single BS)
    par.U = par.MT; % number of single-antenna UEs

%% sim parameters (please read!)
    par.SNRdB_list = [-15:2:25]; % list of SNR [dB] values to be simulated
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
    
%% Testing no rng within channel_sim
% par.UE_x_locs = [52.1595069011440;37.6729557766093;85.7182629858770;36.0725173060400;59.1413893269806;30.0582962630790;40.3291483356594;68.7952297657816;44.2369035243837;79.1779720483239;35.1037738470982;80.6161831832386;42.0630544874686;58.3233107044553;54.9274961614981;85.1932804053457];
% par.UE_y_locs = [5.38010491531581;52.2996911124930;-7.52942946033120;36.1665992379193;11.7439498445019;-43.0907173092592;-51.1208373273156;-46.7836513417716;48.1768049071136;0.502063346574513;-40.3445179555917;-47.1065789042279;-37.9642546794900;8.26654243975507;50.7689197298537;-37.7564019285772];
% par.UE_z_locs = ones(16,1);

par.test = 0;
    
%%
par.iid = 0;            %value of 1: iid Gaussian csi. value of 0: Quadriga
par.array_v = 1;
par.array_h = par.B/par.array_v;
par.shuffle = 0;
decentralized_MIMOsim(par);
hold on
par.shuffle = 1;
decentralized_MIMOsim(par);
% hold on
% par.iid = 0;
% decentralized_MIMOsim(par);
% par.detector = {'DCG'};
% par.iteration = 5;
% par.iid = 0;
% decentralized_MIMOsim(par);
% hold on
% par.iteration = 6;
% 
% decentralized_MIMOsim(par);
% hold off
% legend('uMMSE', 'uMMSE_decent', 'CG', 'DCG_iter4','MF','DCG_iter5', 'DCG_iter6');
% hold on
% par.iid = 0;
% decentralized_MIMOsim(par);
% hold off
% title(['Par.C = ' num2str(par.C) par.scenario ], 'Interpreter','none');