
function par = symbolizer(par)

switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
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
    case '256QAM',
        sym = -15:2:15;
        par.symbols = [];
        for row = 1:16
            for col = 1:16
                par.symbols = [par.symbols sym(row)+sym(col)*1i];
            end
        end
        
    case '1024QAM'
        sym = -31:2:31;
        par.symbols = [];
        for row = 1:length(sym)
            for col = 1:length(sym)
                par.symbols = [par.symbols sym(row)+sym(col)*1i];
            end
        end
    case '4PSK'
        M = 0:3;
        par.symbols = exp(-1i*2*pi*M/4);
    case '8PSK'
        M = 0:7;
        par.symbols = exp(-1i*2*pi*M/8);
    case '32PSK'
        M = 0:31;
        par.symbols = exp(-1i*2*pi*M/32);
    case '16PSK'
        M = 0:15;
        par.symbols = exp(-1i*2*pi*M/16);
    case '64PSK'
        M = 0:63;
        par.symbols = exp(-1i*2*pi*M/64);
    case '256PSK'
        M = 0:255;
        par.symbols = exp(-1i*2*pi*M/256);
    case '1024PSK'
        M = 0:1023;
        par.symbols = exp(-1i*2*pi*M/1024);
    otherwise
        par.symbols = 1;
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);
%% normalize symbol so Es = 1
if par.normalize
    par.symbols = par.symbols/sqrt(par.Es);
    par.Es = 1;
end

par.VarS = par.Es + mean(par.symbols)^2;

end