
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

data_long = readtable('data_long_final.csv');

price_long = table2array(data_long(:,2));
prom_long = table2array(data_long(:,3));
brand_dummies = table2array(data_long(:,5:14));
iv = table2array(data_long(:,15:46));

Z = [iv prom_long brand_dummies];
X1 = [price_long prom_long brand_dummies];

opt_00 = fminsearch(@(params)gmm(params(1:2), X1, Z, mkt_share, prices, income, branded, 50), [0 ; 0]);



%

x = iterate_delta(5, 5, mkt_share, prices, income, branded, 50);







