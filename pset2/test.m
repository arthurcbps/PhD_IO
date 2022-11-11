
data2 = table2array(test);

betas_rho= fminunc(@(params)gmm_crit(params(1:4), Phi, logK, logL, logLagK, logLagL, lagPhi), [0.3; 0.3; 0.4; 0.88]);


logR = log(data2(:, 4));
logK = log(data2(:, 1));
logL = log(data2(:, 2));
logLagK = log(data2(:, 6));
logLagL = log(data2(:, 5));
logLagR = log(data2(:, 8));
lagPhi = data2(:, 10);
Phi = data2(:, 9);


pars = fsolve(@(params)back_out_pars(betas_rho(1:3), params(1), params(2), params(3)), [0.3 0.2 0.8]);


bb= fminunc(@(params)gmm_bb(params(1:4), logR, logK, logL, logLagK, logLagL, logLagR), [0.3; 0.3; 0.4; 0.88]);


gmm_bb(zeros(1,4), logR, logK, logL, logLagK, logLagL, logLagR)

pars_bb = fsolve(@(params)back_out_pars(bb(1:3), params(1), params(2), params(3), params(4)), [0.3 0.2 0.8 -5]);
