function system = back_out_pars(beta, alphal, alphak, theta, eta)
alpham = 1-alphal - alphak;
beta0 = beta(1);
betak = beta(2);
betal = beta(3);

eq1 = beta0 - (1/(1-alpham*theta))*(alpham*theta*log(alpham*theta) - 0.1*theta^2);
eq2 = betal - alphal*theta/(1-alpham*theta);
eq3 = betak - alphak*theta/(1-alpham*theta);
eq4 = theta - (1+eta)/eta;
system = [eq1 eq2 eq3 eq4];

end