
% BLP
clear clc
% load data in wide format (transformed in R script) and extracting
% relevant variables - this is useful for obtaining delta terms

data_wide=readtable('data_clean_wide.csv');


week = table2array(data_wide(:, 1));
store= table2array(data_wide(:, 2));
prices = table2array(data_wide(:, 3:13));
prom = table2array(data_wide(:, 14:24));
income = table2array(data_wide(:, 25:35));
mkt_share = table2array(data_wide(:, 36:46));
branded = table2array(data_wide(:, 47:57));



% load long data for estimation:

data_long = readtable('data_BLPlong.csv');

price_long = table2array(data_long(:,4));
prom_long = table2array(data_long(:,5));
brand_dummies = table2array(data_long(:,6:9));
iv = table2array(data_long(:,10:42));

Z = [iv prom_long brand_dummies];
clear iv
X1 = [price_long prom_long brand_dummies];

% gets nonlinear parameters - sigmas

W = inv(Z'*Z);



[nl_par00, gmm_crit_00] = fminsearch(@(params)gmm(params(1:2), X1, Z, mkt_share, prices, income, branded, 20), [0 ; 0]);
delta_00= iterate_delta(nl_par00(1),nl_par00(2),  mkt_share, prices, income, branded, 50);
delta_00_long= delta_00(:);
clear delta_00
linpar_00 = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_00_long;

%changing the starting point
[nl_par11, gmm_crit_11] = fminsearch(@(params)gmm(params(1:2), X1, Z, mkt_share, prices, income, branded, 20), [1 ; 1]);
delta_11= iterate_delta(nl_par11(1), nl_par11(2),  mkt_share, prices, income, branded, 50);
delta_11_long= delta_11(:);
clear delta_11
linpar_11 = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_11_long;


%% Elasticities for store 9 in week 10

% get data

price_s9_w10 = data_long.price(data_long.store == 9 & data_long.week == 10);
mkt_share_s9_w10 = data_long.mkt_share(data_long.store == 9 & data_long.week == 10);
branded_s9_10 = data_long.branded(data_long.store == 9 & data_long.week == 10);


% first we start with logit, meaning we get deltas fixing sigmas to be 0

% note that it is fine to use the iteration function, logit is nested in
% there
delta_logit = iterate_delta(0, 0, mkt_share, prices, income, branded, 20);
delta_logit_long = delta_logit(:);
% get parts associated with w10 s9

delta_logit_s9w10 = delta_logit(week == 10 & store == 9, :)';

% use the gmm estimator to get alpha

lin_par_logit = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_logit_long;
alpha_logit = lin_par_logit(1);

% now we compute model implied market shares

share = exp(delta_logit_s9w10)./(1 + sum(exp(delta_logit_s9w10)));

slutsky_logit = (1./share).*alpha_logit.*price_s9_w10'.*(diag(share) - share*share');




