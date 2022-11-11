function objective = gmm_crit(par, logR, logK, logL, logLagK, logLagL, lagPhi)

beta0 = par(1);
betak = par(2);
betal = par(3);
rho = par(4);

xi = (logR - beta0 - betak*logK - betal*logL- rho.*(lagPhi - beta0 - betak*logLagK - betal*logLagL));

m1 = xi;
m2 = xi.*logK;
m3 = xi.*logL;
m4 = xi.*lagPhi;
%m5 = xi.*logLagL;

M = [m1 m2 m3 m4];
x=sum(M);

objective = x*x';

end


