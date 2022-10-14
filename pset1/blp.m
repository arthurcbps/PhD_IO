
% BLP
clear clc
% load data in wide format (transformed in R script) and extracting
% relevant variables - this is useful for obtaining delta terms

data_wide=readtable('data_clean_wide.csv');


prices = table2array(data_wide(:, 3:13));
prom = table2array(data_wide(:, 14:24));
income = table2array(data_wide(:, 25:35));
mkt_share = table2array(data_wide(:, 36:46));
branded = table2array(data_wide(:, 47:57));


% load long data for estimation:

data_long = readtable('data_BLPlong.csv');

price_long = table2array(data_long(:,1));
prom_long = table2array(data_long(:,2));
brand_dummies = table2array(data_long(:,3:6));
iv = table2array(data_long(:,7:39));

Z = [iv prom_long brand_dummies];
X1 = [price_long prom_long brand_dummies];

% gets nonlinear parameters - sigmas

W = inv(Z'*Z);



nl_par00 = fminsearch(@(params)gmm(params(1:2), X1, Z, mkt_share, prices, income, branded, 20), [0 ; 0]);
delta_00= iterate_delta(opt_00(1),opt_00(2),  mkt_share, prices, income, branded, 50);
delta_00_long= delta_00(:);

linpar_00 = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_00_long;

%changing the starting point
nl_par11 = fminsearch(@(params)gmm(params(1:2), X1, Z, mkt_share, prices, income, branded, 20), [1 ; 1]);
delta_11= iterate_delta(opt_00(1),opt_00(2),  mkt_share, prices, income, branded, 50);
delta_11_long= delta_11(:);

linpar_11 = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_11_long;











