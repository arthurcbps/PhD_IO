function objective = gmm_bb(par, logR, logK, logL, logLagK, logLagL, logLagR)

beta0 = par(1);
betak = par(2);
betal = par(3);
rho = par(4);

xi = logR - rho*logLagR - beta0*(1-rho) - betak*(logK - rho*logLagK) - betal*(logL - rho*logLagL);

m1 = xi;
m2 = xi.*logK;
m3 = xi.*logL;
m4 = xi.*logLagK;
m5 = xi.*logLagL;
%m5 = xi.*logLagL;

M = [m1 m2 m3 m4 m5];
x=sum(M);

objective = x*x';

end