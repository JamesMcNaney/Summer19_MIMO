% =========================================================================
% -- Simulator for Decentralized Feedforward Massive MU-MIMO Precoding  
% -------------------------------------------------------------------------
% -- (c) 2018 Christoph Studer; e-mail: studer@cornell.edu 
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our
% -- conference paper:
% --    XXXX
% -- and clearly mention this in your paper
% =========================================================================

function dprecoder_sim(varargin)

% add paths that are required
addpath('matlab2tikz');
addpath('param')

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.runId = 0;              % simulation ID (used to reproduce results)
    par.U = 8;                 % number of single-antenna users
    par.B = 128;                % number of base-station antennas (B>>U)
    par.T = 10;                  % number of time slots
    par.C = 2;                  % number of clusters
    par.S = par.B/par.C;
    par.mod = '16QAM';          % modulation type: 'BPSK','QPSK','16QAM','64QAM','8PSK'
    par.trials = 1e2;           % number of Monte-Carlo trials (transmissions)
    par.NTPdB_list = -16:2:14;  % list of normalized transmit power [dB] values
    par.rho2 = 1;               % rho^2=1 (should NOT affect your results!)
    %par.precoder = {'MRT','SMRT','ZF','WF','PD_WF','FD_WF','DP_legacy'};    
    par.precoder = {'MRT','WF','FD_WF'}; 
    par.channel = 'quadriga';   % channel model 'rayleigh', 'los', 'cellfree' 'quadriga'
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
    par.fc = 16e6; % carrier frequency [Hz]
    par.BW = 10e6; % bandwidth [Hz]
    par.N = 1024; % number of carriers
%     par.B = par.MR; % number of antennas in the BS (we use a single BS)
%     par.U = par.MT; % number of single-antenna UEs
    
    
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1};   % only argument is par structure
    
end

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runId);



% simulation name (used for saving results)
par.simName = ['ERR_',num2str(par.U),'x',num2str(par.B), '_C', ...
    num2str(par.C), '_', par.betaest, '_', par.mod, '_', num2str(par.trials),'Trials'];

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    case '8PSK',
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ];
    case '16PSK',
        par.symbols = [ ...
            exp(1i*2*pi*0/16)  ... % 0000
            exp(1i*2*pi*1/16)  ... % 0001     
            exp(1i*2*pi*3/16)  ... % 0010
            exp(1i*2*pi*2/16)  ... % 0011
            exp(1i*2*pi*7/16)  ... % 0100
            exp(1i*2*pi*6/16)  ... % 0101
            exp(1i*2*pi*4/16)  ... % 0110
            exp(1i*2*pi*5/16)  ... % 0111
            exp(1i*2*pi*15/16) ... % 1000
            exp(1i*2*pi*14/16) ... % 1001
            exp(1i*2*pi*12/16) ... % 1010
            exp(1i*2*pi*13/16) ... % 1011
            exp(1i*2*pi*8/16)  ... % 1100
            exp(1i*2*pi*9/16)  ... % 1101
            exp(1i*2*pi*11/16) ... % 1110
            exp(1i*2*pi*10/16) ];  % 1111           
end

% compute symbol energy
par.Es = mean(abs(par.symbols).^2);

% number of antenns per cluster
par.S = par.B/par.C;

% precompute bit labels
par.bps = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.bps,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% - initialize result arrays (detector x normalized transmit power)
[res.PER, res.SER, res.BER ] = deal(zeros(length(par.precoder),length(par.NTPdB_list)));
[res.TxPower, res.RxPower, res.TIME] = deal(zeros(length(par.precoder),length(par.NTPdB_list)));

% compute noise variances to be considered: NTP = rho^2/N0
N0_list = par.rho2*10.^(-par.NTPdB_list/10);


% trials loop
tic
for t=1:par.trials
    
    % generate data
    for qq=1:par.T
        % generate random bit stream
        B(:,:,qq) = randi([0 1],par.U,par.bps);
        % generate transmit symbol
        Idx(:,qq) = bi2de(B(:,:,qq),'left-msb')+1;
        S(:,qq) = par.symbols(Idx(:,qq)).';
    end
    
    % generate masks
    MaskI = true(par.U,par.T);
    MaskT = true(1,par.T);
    switch par.betaest
        case {'pilot'}
            S(:,1) = ones(par.U,1)*sqrt(par.Es); % just send all ones in the first time slots
            MaskI(:,1) = false;
            MaskT(1,1) = false;
    end
    
    
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
    
    % algorithm loop
    for d=1:length(par.precoder)
        
        % normalized transmit power loop
        for k=1:length(par.NTPdB_list)
            
            % set noise variance
            N0 = N0_list(k);
            
            % record time used by the beamformer
            starttime = toc;
            
            % beamformers
            switch (par.precoder{d})
                case 'MRT',     % MRT beamforming
                    [X, beta] = MRT(par,S,H,N0);     
                case 'SMRT',     % scaled MRT beamforming
                    [X, beta] = SMRT(par,S,H,N0);                      
                case 'ZF',      % ZF beamforming
                    [X, beta] = ZF(par,S,H,N0);              
                case 'CD',      % CD beamforming
                    [X, beta] = CD(par,S,H,N0);
                case 'DCD',      % ZF beamforming
                    [X, beta] = DCD(par,S,H,N0); 
                case 'WF',
                    [X, beta] = WF(par,S,H,N0);
                case 'PD_WF',
                    [X, beta] = PD_WF(par,S,H,N0);                        
                case 'FD_WF',
                    [X, beta] = FD_WF(par,S,H,N0);         
                case 'DP_legacy',
                    [X, beta] = DP_legacy(par,S,H,N0);                         
                otherwise,
                    error('par.precoder not specified')
            end
            
            % record beamforming simulation time
            res.TIME(d,k) = res.TIME(d,k) + (toc-starttime);
            
            
            % transmit data over noisy channel
            HX = H*X;
            Y = HX + sqrt(N0)*N;
            
            % extract transmit and receive power
            res.TxPower(d,k) = res.TxPower(d,k) + mean(sum(abs(X(:)).^2))/par.T;
            res.RxPower(d,k) = res.RxPower(d,k) + mean(sum(abs(HX(:)).^2))/par.U/par.T;
            
            % UEs must estimate beta
            switch par.betaest
                case 'genie', % perfect beta directly from beamformer
                    betaest = ones(par.U,1)*beta;
                case 'pilot', % knows that first symbols are for training
                    betaest = real(1./Y(:,1)*sqrt(par.Es)); % ML estimate since we have no prior on beta
            end
            
            % perform estiamtion
            Shat = (betaest*ones(1,par.T)).*Y;
            
            % UE-side detection
            for qq=1:par.T
                [~,Idxhat(:,qq)] = min(abs(Shat(:,qq)*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
                Bhat(:,:,qq) = par.bits(Idxhat(:,qq),:);
            end
            
            % -- compute error and complexity metrics
            err = (Idx(MaskI)~=Idxhat(MaskI));
            res.PER(d,k) = res.PER(d,k) + any(err(:));
            res.SER(d,k) = res.SER(d,k) + sum(err(:))/par.U/par.T;
            tmpBER = B(:,:,MaskT)~=Bhat(:,:,MaskT);
            res.BER(d,k) = res.BER(d,k) + sum(tmpBER(:))/(par.U*par.bps*sum(MaskT));
                        
        end % NTP loop
        
    end % algorithm loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',...
            time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.PER = res.PER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.TxPower = res.TxPower/par.trials;
res.RxPower = res.RxPower/par.trials;
res.TIME = res.TIME/par.trials;


%res.TIME 
%res.RxPower
res.TxPower

% -- save final results (par and res structures)

if par.save
    save([ 'results/' par.simName '_' num2str(par.runId) ],'par','res');
end

% -- show results (generates fairly nice Matlab plots)

if par.plot
    
    % - BER results
    
  if(par.channel == 'quadriga')
    marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  else
    marker_style = {'go-','cs--','kv-.','kp:','g*-','c>--','yx:'};
  end
    
%     marker_style = {'kx-','bo:','rs--','mv-.','gp-.','bs--','y*--'};
    h = figure(1);
    for d=1:length(par.precoder)
        semilogy(par.NTPdB_list,res.BER(d,:),marker_style{d},'LineWidth',2);
        if (d==1)
            hold on
        end
    end
    hold on
    grid on
    box on
    xlabel('normalized transmit power [dB]','FontSize',12)
    ylabel('uncoded bit error rate (BER)','FontSize',12);
    if length(par.NTPdB_list) > 1
        axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-3 1]);
    end
    if (par.iid == 0)
        legend([strcat('rayleigh', '_', par.precoder) strcat(par.channel, '_', par.precoder)],'FontSize',12,'location','southwest','Interpreter','none')
        set(gca,'FontSize',12);
    end
    if par.save
        % save eps figure (in color and with a reasonable bounding box)
        print(h,'-loose','-depsc',[ 'results/',par.simName '_' num2str(par.runId) ])
        matlab2tikz(['results/',par.simName,'_' num2str(par.runId), '.tex'],'standalone',true);
    end
    
end

end


%% Maximum ratio transmission (MRT) beamforming
function [X, beta, P] = MRT(par,S,H,N0)

% transmitted signal
P = H';
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));

%G = H*H';
%betainv = sqrt(par.rho2)/sqrt(par.Es*sum(1./diag(G)));

X = betainv*(P*S);

% average scaling over signals
beta = 1/betainv;

%beta = 1.0*(norm(s,2)^2+N0*par.U)/(s'*H*x); % that's cheating

end

%% Scaled maximum ratio transmission (SMRT) beamforming
function [X, beta, P] = SMRT(par,S,H,N0)

% transmitted signal
G = H*H';
P = H'*diag(1./diag(G));

%G = H*H';
betainv = sqrt(par.rho2)/sqrt(par.Es*sum(1./diag(G)));
%betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));

X = betainv*(P*S);

%X = P*S;
%betainv = sqrt(par.rho2)/norm(X,2);
%X = betainv*X;

% average scaling over signals
beta = 1/betainv;

end


%% Zero-forcing (ZF) beamforming
function [X, beta] = ZF(par,S,H,N0)

% transmitted signal
P = zfinv(par,H);
%betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));

G = H*H';
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(inv(G)));

X = betainv*(P*S);

% average scaling over signals
beta = 1/betainv;

end

% ZF pseudo inverse
function Hinv = zfinv(par,H)
[U,S] = size(H);
if S>=U
    Hinv = H'/(H*H');
else
    Hinv = (H'*H)\H';
end
end

function [X, beta] = CD(par,S,H,N0)

[row, col]=size(H);

% instant normalization constant
%rho = sqrt((par.B-par.U)/(par.Es*par.U));
%rho = 1/sqrt(par.Es*trace(inv(H*H')));

dinv = zeros(row, 1);
H_norm = zeros(row, col); 
s_norm = zeros(row);

X = zeros(col, par.T);

for t=1:par.T
    % -- preprocessing
    for u=1:row
        dinv(u) = 1/norm(H(u,:),2);
        H_norm(u, :) = dinv(u)*H(u,:);
        s_norm(u) = dinv(u)*S(u, t);
    end
    % -- CD iterations
    for k=1:par.maxiter
        for u=1:row
           X(:, t) = X(:, t) - par.steplen*(H_norm(u,:) * X(:, t) - s_norm(u))* H_norm(u, :)';
        end
    end
end

%A = inv(H*H');
%betainv = sqrt(par.rho2)/(sqrt(par.Es*trace(A)));

%G=H*H';
%betainv = sqrt(par.rho2)/(sqrt(par.Es*sum(1./eig(G))));
%G=H_norm*H_norm';
%betainv = sqrt(par.rho2)/(sqrt(par.Es*sum(eig(G))));


betainv = par.damping * sqrt(par.rho2)/(sqrt(par.Es*sum(dinv.^2)));
X =  betainv*X;

%X = sqrt(par.rho2)*X/norm(X,2);
%betainv=sqrt(par.rho2)/norm(X,2);

beta = 1/betainv;


end

function [xout, betainv] = DCD_solver(par, s, H)

[row, col]=size(H);

% -- preprocessing
dinv = zeros(row, 1);
x = zeros(col, 1);

H_norm = zeros(row, col); 
s_norm = zeros(row);

for u=1:row
    dinv(u) = 1/norm(H(u,:),2);
    H_norm(u, :) = dinv(u)*H(u,:);
    s_norm(u) = dinv(u)*s(u);
end

% -- CD iterations
for k=1:par.maxiter
    for u=1:row
       x = x - par.steplen*(H_norm(u,:) * x - s_norm(u))* H_norm(u, :)';
    end
end

%A = inv(H*H');
%betainv = sqrt(par.rho2/par.C)/(sqrt(par.Es*trace(A)));
betainv = par.damping*sqrt(par.rho2/par.C)/(sqrt(par.Es*sum(dinv.^2)));
    
xout = betainv*x;

end

function [X, beta] = DCD(par,S,H,N0)

% extract per-cluster channel matrix
for cc=1:par.C
    Hc(:,:,cc) = H(:,(cc-1)*par.S+1:cc*par.S);     
end

% beamformer
betainv=zeros(1, par.C);
xc = zeros(par.S, par.C, par.T);
for t=1:par.T
    for cc=1:par.C
        [xc(:, cc, t), betainv(:, cc)] = DCD_solver(par, S(:, t), Hc(:,:,cc));
    end
end

X = reshape(xc, par.S*par.C, par.T);

% transmitted signal

beta = 1/sum(betainv);

end


%% Centralized Wiener-filter (WF) beamforming
function [X, beta] = WF(par,S,H,N0)

% compute regularized inverse
Ainv = rpinv(par,H,par.U*N0/par.rho2); 

% calculate precoding factor efficiently
beta = sqrt((par.Es/par.rho2)*(trace(Ainv)-sum(abs(Ainv(:)).^2)*(par.U*N0/par.rho2)));

% apply inverse to data in centralized manner
if par.B>=par.U
    X = (1/beta)*(H'*(Ainv*S));
else
    X = (1/beta)*(Ainv*(H'*S));
end

end


% Wiener filter regularized pseudo inverse
function Ainv = rpinv(par,H,reg)
[U,S] = size(H);
if S>=U
    Ainv = inv(H*H'+(reg)*eye(par.U));
else
    Ainv = inv(H'*H+(reg)*eye(par.S));
end
end


%% Partially-Decentralized Wiener-filter (WF) beamforming
function [X, beta] = PD_WF(par,S,H,N0)

% decentralized gram matrix computation and localized averaging
Gc = zeros(par.U,par.U);
for cc=1:par.C
    % calculate local Gram matrix
    Hc(:,:,cc) = H(:,(cc-1)*par.S+1:cc*par.S);
    % average among clusters (can be done in a tree-like fashion)
    Gc = Gc + Hc(:,:,cc)*Hc(:,:,cc)'; 
end

% compute whitening filter at centralized node
Ainv = inv(Gc+(par.U*N0/par.rho2)*eye(par.U));

% calculate precoding factor efficiently
beta = sqrt((par.Es/par.rho2)*(trace(Ainv)-sum(abs(Ainv(:)).^2)*(par.U*N0/par.rho2)));

% whiten transmit signals
Z = (1/beta)*(Ainv*S);

% perform decentralized MRT with whitened signals
for cc=1:par.C
    X((cc-1)*par.S+1:cc*par.S,:) = Hc(:,:,cc)'*Z;    
end

end

%% Fully-Decentralized Wiener-filter (WF) beamforming
function [X, beta] = FD_WF(par,S,H,N0)

%initialization 
stomp = par.FD_WF.stomp;

% perform fully decentralized WF precoding 
for cc=1:par.C
    
    % calculate local precoding matrix
    Hc = H(:,(cc-1)*par.S+1:cc*par.S);
    Ainvc = rpinv(par,Hc,stomp*par.U*N0/par.rho2);    
    betac(cc,1) = sqrt((par.Es/(par.rho2/par.C))*(trace(Ainvc)-sum(abs(Ainvc(:)).^2)*(stomp*par.U*N0/(par.rho2))));
    
    %betainv(cc,1) = sqrt(par.rho2/par.C)/sqrt(par.Es*trace(Pc*Pc'));
    
    if par.S>=par.U    
        X((cc-1)*par.S+1:cc*par.S,:) = (1/betac(cc,1))*(Hc'*(Ainvc*S));
    else
        X((cc-1)*par.S+1:cc*par.S,:) = (1/betac(cc,1))*(Ainvc*(Hc'*S));        
    end
    
end

% average scaling over signals
beta = sum(betac);
%beta = 1/(sum(1./betac));

end



%% decentralized precoder, ADMM version
function [X,beta] = DP_legacy(par,S,H,N0)

%  -- initialize
delta = par.DP_legacy.delta;
gamma = par.DP_legacy.gamma;
maxiter = par.DP_legacy.maxiter;

H_c = zeros(par.U,par.S,par.C);
AinvH = zeros(par.S,par.U,par.C);

% -- preprocessing
for c=1:par.C
    H_c(:,:,c) = H(:,par.S*(c-1)+1:par.S*c); % get the appropriate part of H
    AinvH(:,:,c) = (H_c(:,:,c)'*H_c(:,:,c) + (1/delta)*(par.U*N0/par.rho2)*eye(par.S))\H_c(:,:,c)';
end

for tt=1:par.T

    % initialize running variables
    lambda_c = zeros(par.U,par.C);
    x_c = zeros(par.S,par.C);
    w_c = zeros(par.U,par.C);
    Hx_c = zeros(par.U,par.C);
        
    s = S(:,tt);
    
    % important for fast convergence (reasonable initial guess)
    z_c = max(par.U/par.B,1/par.C)*s*ones(1,par.C);    
    
    % -- start iteration
    for ll = 1:maxiter
    
        % cluster-wise equalization
        for c=1:par.C
            x_c(:,c) = AinvH(:,:,c)*(z_c(:,c) + lambda_c(:,c)); 
            Hx_c(:,c) = H_c(:,:,c)*x_c(:,c);
            w_c(:,c) = (Hx_c(:,c)-lambda_c(:,c));
        end
                   
        % consensus step
        w_avg = 1/(par.C*delta+delta^2)*(par.C*s+delta*sum(w_c,2));
    
        % cluster-wise update
        for c=1:par.C
            z_c(:,c) = (1/delta)*(s+delta*w_c(:,c))-w_avg;
            lambda_c(:,c) = lambda_c(:,c) - gamma*(Hx_c(:,c)-z_c(:,c));
        end
    end
    
    X(:,tt) = x_c(:); % vectorize output
    
end

% instantaneous power normalization
X = X*(sqrt(par.rho2)/sqrt(sum(abs(X(:)).^2)/par.T));

beta = NaN;

end
