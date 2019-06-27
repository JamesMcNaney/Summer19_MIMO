function [res] = detector(res,par,H,y,N0,s,idx,k,t,bits)
% algorithm loop
for d=1:length(par.detector)
    
    switch (par.detector{d}) % select algorithms
        case 'MF'
            [idxhat,bithat] = MF(par,y,H);
        case 'SIMO'
            [idxhat,bithat] = SIMO_bound(par,y,H,s);
        case 'AWGN'
            [idxhat,bithat] = SISO_AWGN(par,N0,s);
        case 'ZF', % zero-forcing detection
            [idxhat,bithat] = ZF(par,H,y,s);
        case 'bMMSE', % biased MMSE detector
            [idxhat,bithat] = bMMSE(par,H,y,N0);
        case 'uMMSE', % unbiased MMSE detector
            [idxhat,bithat] = uMMSE(par,H,y,N0,s);
        case 'uMMSE_decent', % unbiased MMSE detector
            [idxhat,bithat] = uMMSE_decent(par,H,y,N0);
        case 'CG', % CG detector
            [idxhat,bithat] = CG(par,H,y,N0);
        case 'DCG', % decentralized CG detector
            [idxhat,bithat] = DCG(par,H,y,N0);
        otherwise
            disp('no detectors!')
    end
    
    
    % -- compute error metrics
    err = (idx~=idxhat);
    res.VER(d,k) = res.VER(d,k) + any(err);
    res.SER(d,k) = res.SER(d,k) + sum(err)/par.MT;
    res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);
    res.ERR(d,k) = res.ERR(d,k) + norm(y-H*(par.symbols(idxhat).'))^2/par.MR;
    
    if strcmp(par.detector{d},'IO-LAMA')
        for iter = 1:par.iteration
            %% IO
            err = (idx~=idxhat_vec(:,iter));
            res.LAMA_iter_VER(iter,k) = res.LAMA_iter_VER(iter,k) + any(err);
            res.LAMA_iter_SER(iter,k) = res.LAMA_iter_SER(iter,k) + sum(err)/par.MT;
            res.LAMA_iter_BER(iter,k) = res.LAMA_iter_BER(iter,k) + sum(sum(bits(:,:,t)~=bithat_vec(:,:,iter)))/(par.MT*par.Q);
        end
    end
    
end % algorithm loop

end

% -- detector functions

% MF
function [idxhat,bithat] = MF(par,y,H)
xhat = H'*y;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

% SIMO bound
function [idxhat,bithat] = SIMO_bound(par,y,H,s)
n = y - H*s;
idxhat = zeros(par.MT,1);
for idx=1:par.MT
    y_i = H(:,idx)*s(idx) +n;
    y_MF = H(:,idx)'*y_i/norm(H(:,idx),2)^2;
    [~, idxhat(idx)] = min(abs(y_MF*ones(1,length(par.symbols))-par.symbols).^2,[],2);
end
bithat = par.bits(idxhat,:);
end

% SISO channel detector
function [idxhat,bithat] = SISO_AWGN(par,N0,s)
xhat = s+sqrt(0.5*(N0))*(randn(par.MT,1)+1i*randn(par.MT,1));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end


% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y,sN)
xhat = H\y;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%MMSE
function [idxhat,bithat] = MMSE(par,H,y,N0)
W = (H'*H+(N0/par.Es)*eye(par.MT))\(H');
xhat = W*y;
G = real(diag(W*H));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

% unbiased MMSE detector (uMMSE)
function [idxhat,bithat,xhat] = uMMSE(par,H,y,N0,sN)

%% tse-hanly
%
% b = -N0 + par.Es*(1-par.beta);
% sigsq_hat_cen = (-b + ...
%     sqrt(b.^2 + 4*N0*par.Es))/2
% var(xhat-sN)
%
W = (H'*H+(0.5*N0/par.Es)*eye(par.MT))\(H');
xhat = W*y;
G = real(diag(W*H));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end


% unbiased MMSE detector (uMMSE)
function [idxhat,bithat,xhat] = uMMSE_decent(par,H,y,N0)

H_decent = zeros(par.MR/par.C,par.MT,par.C);
idx_vec = reshape(1:par.MR,par.MR/par.C,par.C);

for m=1:par.C
    H_decent(:,:,m) = H(idx_vec(:,m)',:);
end



y_decent = reshape(y,par.MR/par.C,par.C);

xhat = zeros(par.MT,par.C);
for j=1:par.C
    W = (H_decent(:,:,j)'*H_decent(:,:,j)+(0.5*N0/par.Es)*eye(par.MT))\(H_decent(:,:,j)');
    xhat(:,j) = W*y_decent(:,j);
end
xhat_final = mean(xhat.').';

W = (H'*H+(0.5*N0/par.Es)*eye(par.MT))\(H');
G = real(diag(W*H));
[~,idxhat] = min(abs(xhat_final*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end


function [idxhat,bithat] = CG(par,H,y,N0)

r = CGDET(par, H, y, N0);
[~, idxhat]=min(abs(r*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end

function [idxhat,bithat] = DCG(par,H,y,N0)

% decent_channel
H_decent = zeros(par.MR/par.C,par.MT,par.C);
idx_vec = reshape(1:par.MR,par.MR/par.C,par.C);

for m=1:par.C
    H_decent(:,:,m) = H(idx_vec(:,m)',:);
end

y_decent = reshape(y,par.MR/par.C,par.C);

r = zeros(par.MT,par.C);

par_decent = par;
par_decent.MR = par.MR/par.C;
par_decent.beta = par.beta*par.C;

for c=1:par.C
    r(:,c) = CGDET(par_decent,H_decent(:,:,c),y_decent(:,c),N0/par.C);   
end

% average
r_decent = sum(r,2)/par.C;

[~, idxhat]=min(abs(r_decent*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);

%W = (H'*H+(N0/par.Es)*eye(par.MT))\(H');
%G = real(diag(W*H));
%[~,idxhat] = min(abs(r_decent*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);

bithat = par.bits(idxhat,:);
end

