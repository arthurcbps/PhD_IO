
% BLP
clear clc
% load data in long formatted (transformed in R script) and extracting
% relevant variables

data_wide=readtable('data_clean_wide.csv');


prices = table2array(data_wide(:, 3:13));
prom = table2array(data_wide(:, 14:24));
income = table2array(data_wide(:, 25:35));
mkt_share = table2array(data_wide(:, 36:46));
branded = table2array(data_wide(:, 47:57));

clear data_long


% first step - get implied model shares
% draws
R = 50;

x = iterate_delta(5, 5, mkt_share, prices, income, branded);







