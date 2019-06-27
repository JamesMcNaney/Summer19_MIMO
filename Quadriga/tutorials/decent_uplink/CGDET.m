function  res = CGDET(par,H,y,N0)
% -- preprocessing
A = (H'*H+(0.5*N0/par.Es)*eye(par.MT)); % MMSE

b=H'*y;
r0 = b;
p0 = r0;

rk1 = r0;
xk1 = zeros(par.MT,1);
rk1_sq = r0'*r0;
pk1 = p0;

alphak = zeros(1,par.iteration);
betak = zeros(1,par.iteration);

for k = 1:par.iteration
    Apk1 = A*pk1;
    alphak(k) = rk1_sq/(pk1'*Apk1);
    
    alphak_pk1 = alphak(k)*pk1;
    xk = xk1 + alphak_pk1;
    
    rk = rk1 - alphak(k)*Apk1;
    rk_sq = rk'*rk;
    betak(k) = (rk_sq/rk1_sq);
    pk = rk + betak(k)*pk1;
    
    rk1_sq = rk_sq;
    xk1= xk;
    pk1 = pk;
    rk1 = rk;    
end

res = xk;

end