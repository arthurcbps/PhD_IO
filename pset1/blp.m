
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

x = iterate_delta(1, 1, mkt_share, prices, income, branded);




delta0 = zeros(size(prices));
delta_hat = log(mkt_share) - log(share_fun(delta0, income, branded, prices));

tol = 1e-6;
norm = 1;

i=0;
while norm > tol
    delold = delta_hat;
    delta_hat = delold + log(mkt_share) - log(share_fun(delold, income, branded, prices));
    
    t = abs(delta_hat - delold);
    norm = max(t, [], 'all');
    i = i + 1;
disp(['error:' num2str(norm)])

end

x= log(share_fun(delold, income, branded, prices));


