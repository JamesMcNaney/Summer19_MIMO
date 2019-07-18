clear all;

par.runId = 0;              % simulation ID (used to reproduce results)
par.U = 8;                 % number of single-antenna users
par.B = 128;                % number of base-station antennas (B>>U)
par.T = 10;                  % number of time slots
par.C = 8;                  % number of clusters
par.S = par.B/par.C;
par.mod = '16QAM';          % modulation type: 'BPSK','QPSK','16QAM','64QAM','8PSK'
par.trials = 2e2;           % number of Monte-Carlo trials (transmissions)
par.NTPdB_list = -16:2:14;  % list of normalized transmit power [dB] values
par.rho2 = 1;               % rho^2=1 (should NOT affect your results!)
%par.precoder = {'MRT','SMRT','ZF','WF','PD_WF','FD_WF','DP_legacy'};    
par.precoder = {'PD_WF'}; 
par.channel = 'quadriga';   % channel model 'rayleigh', 'los', 'cellfree' 'quadriga'
par.iid = 1;
par.betaest = 'pilot';      % 'pilot', 'genie'
par.save = false;           % save results (true,false)
par.plot = true;            % plot results (true,false)
% algorithm-dependent parameters
% FD_WF
par.FD_WF.stomp = 0.125*par.C; % determines how much to stomp regularization (C)
% DP_legacy
par.DP_legacy.delta = 0.3; % Lagrange scaling (1.0)
par.DP_legacy.gamma = 1.0; % Lagrange stepsize (1.0)
par.DP_legacy.maxiter = 2; % keep at 2

par.maxiter = 4; % max number of iterations for OCD algorithms
par.steplen = 1; % step size for CD
par.damping = 1; % damping CD energy caused by approximation

%to be deleted/changed, trying to integrated quadriga
par.scenario = 'Freespace'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
par.fc = 3.5e9; % carrier frequency [Hz]
par.BW = 10e6; % bandwidth [Hz]
par.N = 1024; % number of carriers


par.array_v = 1;        % number of vertical antenna elements, assuming total
                    % number of antenna elements is a power of 2, must also
                    % be a power of 2 < total number of antenna elements

par.array_h = par.B/par.array_v;

par.iter = 0;
downlink(par);
hold on
par.array_v = 2;
par.array_h = par.B/par.array_v;
par.iter = 1;
downlink(par);
hold on
par.array_v = 4;
par.array_h = par.B/par.array_v;
par.iter = 2;
downlink(par);
hold on
par.array_v = 8;
par.array_h = par.B/par.array_v;
par.iter = 3;
par.legend = ["1x128", "2x64", "4x32", "8x16"];
par.iid = 0;
downlink(par);
hold off

% title(strcat(par.scenario, ' users = ', num2str(par.U),' par.mod = ',par.mod, ' par.C =', num2str(par.C)),'Interpreter', 'none');
title(strcat(par.scenario, ' precoder = ', par.precoder, ' par.C =', num2str(par.C)),'Interpreter', 'none');