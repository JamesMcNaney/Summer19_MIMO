function decentralized_MIMOsim(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    clc
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.runId = 0; % simulation ID (used to reproduce results)
    par.C = 2;
    par.MR = 128; % total receive antennas
    
    par.MT = 16; % transmit antennas (set not larger than MR!)
    
    par.beta = par.MT/par.MR;
    % number of clusters (uniformly distributed antennas over C)
    % (set C>1 as C=1 gives an error for the decentralized detectors
    % MR/C has to be integer!

    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 5000; % number of Monte-Carlo trials (transmissions)    
    
%% added parameters for QuaDRiGa
    par.scenario = 'LOSonly'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
    par.fc = 16e6; % carrier frequency [Hz]
    par.BW = 10e6; % bandwidth [Hz]
    par.N = 1024; % number of carriers
    par.B = par.MR; % number of antennas in the BS (we use a single BS)
    par.U = par.MT; % number of single-antenna UEs
%     par.iid = 1;
%% sim parameters (please read!)
    par.SNRdB_list = [-15:5:40]; % list of SNR [dB] values to be simulated
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
else
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % argument is par structure
end

% -- initialization

% use runId random seed (improved reproducibility)
rng(par.runId,'twister');

% Define symbols
par = symbolizer(par);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

% -- start simulation

% initialize result arrays (detector x SNR x EVM)
res.VER = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
res.SER = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate

% error
res.ERR = zeros(length(par.detector),length(par.SNRdB_list));

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.MT,par.Q,par.trials);

% loads in the Quadriga channel batch of coefficient matrices
% H_batch = load(strcat(par.scenario, '_linearAntenna_128B_16U_4000batch.mat'));
if par.array_v == 4
    H_batch = load(strcat(par.scenario, '_4x32Antenna_128B_16U_4000batch.mat'));
else
    H_batch = load(strcat(par.scenario, '_linearAntenna_128B_16U_4000batch.mat'));
end
brute_force_shuffle_idx = [1:8,33:40,65:72,97:104,9:16,41:48,73:80,105:112,17:24,49:56,81:88,113:120,25:32,57:64,89:96,121:128]';

% trials loop
for t=1:par.trials
    complete = 100*t/par.trials;
    disp([num2str(complete,'%.3f') '% completed'])
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    %         n = randn(par.MR,1);
    % LAMA channel matrix - gaussian with variance 1/MR
    if par.iid == 1
      H = sqrt(0.5/par.MR)*...
          (randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    else  
        H = H_batch.csi_batch(:,:,(randi([1,4000])));
%       H = channel_sim(par);
%       H = normalize(H);
%         norm_coef = zeros(1,par.MT);                    
%         for i = 1:par.MT
%             for j = 1:par.MR
%                 norm_coef(i)=norm_coef(i)+norm(H(j,i)); %sum the 2-norms of each column
%             end
%             norm_coef(i) = norm_coef(i)/par.MR;         %average the 2-norm sum
%             H(:,i) = H(:,i)/norm_coef(i)*.07833; 
%         end
%         accum  = zeros(par.MR,par.MT);
%       for i = 1:par.MR
%           for j = 1:par.MT
%               accum(i,j) = norm(H(i,j));
%           end
%       end
%       accum2 = zeros(1,par.MT);
%       for i = 1:par.MT
%           accum2(i) = sum(accum(:,i))/par.MR;
%       end
%       par.test = par.test + sum(accum2)/par.MT;
%       if t == 99
%           dummy = 1;
%       end
        if par.shuffle == 1
            indexer = reshape(1:par.MR, par.C,par.MR/par.C);
            Hidxer = [];
            for i = 1:par.C
                Hidxer = [Hidxer indexer(i,:)];
            end
            H = H(Hidxer,:);
        end
        if par.shuffle == 2
            H = H(brute_force_shuffle_idx,:);
        end
    end
    % transmit over noiseless channel (will be used later)
    x_send = H*s;
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
        N0 = par.beta*par.Es*10^(-par.SNRdB_list(k)/10);
        
        % transmit data over noisy channel
        y = x_send+sqrt(N0)*n;
        
        % assume pilot tone = identity so, H_est = H + N
        H_N0 = sqrt(0.5*N0)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
        H_est = H + H_N0;
        
        
        %% detection
        [res] = detector(res,par,H,y,N0,s,idx,k,t,bits);
        
    end % SNR loop
end % trials loop


% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.ERR = res.ERR/par.trials;

save([ par.simName '_' num2str(par.runId) ],'par','res');

% -- plot results
% if par.iid == 0
%     marker_style = {'ko-','ro-','bs-','bo-','bv-'};
% else
%     marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:','ko-',...
%     'b^--','r*-.','mo:','kv-','gp--','c*-.','y>:','kx-'};
% end

% if par.iteration == 5
%     marker_style = {'b+-'};
% end
% if par.iteration == 6
%         marker_style = {'bx-'};
% end

% if par.iteration == 4
if par.shuffle == 0
% if(par.iid == 1)
%     marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
%   else
    marker_style = {'go-','co--','ro-.','bo-','ko--','kp:','g*-','c>--','yx:'};
%   end
end

if par.shuffle == 1
    marker_style = {'gx-','cx--','rx-.','bx-','kx--','kp:','g*-','c>--','yx:'};
end

if par.shuffle == 2
    marker_style = {'g*-','c*--','r*-.','b*-','kx--','kp:','g*-','c>--','yx:'};
end

%% figure
for d=1:length(par.detector)
    if d==1
        semilogy(par.SNRdB_list,res.BER(d,:),...
            marker_style{d},'LineWidth',2)
        %                 semilogy(par.SNRdB_list,res.ERR(d,:),...
        %             marker_style{d},'LineWidth',2)
        hold on
    else
        semilogy(par.SNRdB_list,res.BER(d,:),...
            marker_style{d},'LineWidth',2)
        %                 semilogy(par.SNRdB_list,res.ERR(d,:),...
        %             marker_style{d},'LineWidth',2)
    end
end
%hold off
hold on
grid on
xlabel('average SNR per receive antenna [dB]','FontSize',12);
%ylabel('symbol error rate (SER)','FontSize',12);
ylabel('bit error rate (BER)','FontSize',12);
axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-3 1]);
legend(par.detector,'FontSize',12,'location','southwest','interpreter','none');

% if par.iteration == 4
% if(par.iid == 0)
%     legend([strcat('iid_',string(par.detector)),strcat('QuaDRiGa_',string(par.detector))],'Fontsize',12,'Interpreter','none');
%     par.legend = [par.legend string(par.detector)];
%     legend(par.legend,'Fontsize',12,'Interpreter','none');    
%     set(gca,'FontSize',12)
% end
% end

end
