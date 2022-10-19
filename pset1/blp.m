
% BLP
clear clc close all
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
product_dummies = table2array(data_long(:, 43:53));
cost = table2array(data_long(:,10));
price_iv = table2array(data_long(:,12:42));


Z = [cost price_iv prom_long product_dummies];
clear iv
X1 = [price_long prom_long product_dummies];


W = inv(Z'*Z);

% number of draws
R = 20;

%optimization
% gets nonlinear parameters - sigmas

[nl_par00, gmm_crit_00] = fminsearch(@(params)gmm(params(1), params(2), X1, Z, mkt_share, prices, income, branded, R), [0 ; 0]);
delta_00= iterate_delta(nl_par00(1), nl_par00(2), mkt_share, prices, income, branded, R);
delta_00_long= reshape(delta_00', [], 1);
clear delta_00
linpar_00 = inv(X1'*Z*inv(Z'*Z)*Z'*X1)*X1'*Z*inv(Z'*Z)*Z'*delta_00_long;

%changing the starting point
[nl_par11, gmm_crit_11] = fminsearch(@(params)gmm(params(1), params(2), X1, Z, mkt_share, prices, income, branded, R), [1 ; 1]);
delta_11= iterate_delta(nl_par11(1), nl_par11(2),  mkt_share, prices, income, branded, R);
delta_11_long= reshape(delta_11', [], 1);
clear delta_11
linpar_11 = inv(X1'*Z*inv(Z'*Z)*Z'*X1)*X1'*Z*inv(Z'*Z)*Z'*delta_11_long;



% alternative estimation (just in case) - using product dummies instead of
% brand dummies
% iv = table2array(data_long(:,10:42));
% Zalt = [iv prom_long product_dummies];
% clear iv
% X1alt = [price_long prom_long product_dummies];

% gets nonlinear parameters - sigmas

% Walt = inv(Zalt'*Zalt);
% 
% 
% 
% [nl_paralt, gmm_crit_alt] = fminsearch(@(params)gmm(params(1), params(2), X1alt, Zalt, mkt_share, prices, income, branded, R), [0 ; 0]);
% delta_alt= iterate_delta(nl_paralt(1),nl_paralt(2),  mkt_share, prices, income, branded, R);
% delta_alt_long= reshape(delta_alt', [], 1);
% clear delta_alt
% linpar_alt= inv(X1alt'*Zalt*Walt*Zalt'*X1alt)*X1alt'*Zalt*Walt*Zalt'*delta_alt_long;


%% Elasticities for store 9 in week 10

% get data

price_s9_w10 = data_long.price(data_long.store == 9 & data_long.week == 10);
mkt_share_s9_w10 = data_long.mkt_share(data_long.store == 9 & data_long.week == 10);


% first we start with logit, meaning we get deltas fixing sigmas to be 0

% note that it is fine to use the iteration function, logit is nested in
% there
delta_logit = iterate_delta(0, 0, mkt_share, prices, income, branded, R);
delta_logit_long = delta_logit(:);
% get parts associated with w10 s9

delta_logit_s9w10 = delta_logit(week == 10 & store == 9, :)';

% use the gmm estimator to get alpha

lin_par_logit = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_logit_long;
alpha_logit = lin_par_logit(1);

% now we compute model implied market shares (should be the same as actual
% market shares at the true values - so this is just a sanity check)

share = exp(delta_logit_s9w10)./(1 + sum(exp(delta_logit_s9w10)));
% and formula for substitution matrix
slutsky_logit = (1./share).*alpha_logit.*price_s9_w10'.*(diag(share) - share*share');


% Now for BLP elasticities, there are possibly many ways to do so, but I'll
% draw from the income distribution and from the gaussian error. Then we
% compute the Slutsky matrix for each draw and take the average, in the end
% we divide by observed market shares (which of course are not the same as
% model implied ones anymore) - this follow Nevo's approach


% get deltas for week 10, store 9 products
alpha_blp = linpar_00(1);
sigmaI_blp = nl_par00(1);
sigmaB_blp = nl_par00(2);

delta_blp= iterate_delta(sigmaI_blp, sigmaB_blp, mkt_share, prices, income, branded, R);
% get additional parts associated with w10 s9
delta_blp_s9w10 = delta_blp(week == 10 & store == 9, :);
branded_s9w10 = branded(week == 10 & store == 9, :);


rng(1)
v_draws = randn(R,1)';

indices = randi(length(income), 1, R);
income_draws = income(indices, 1);
clear indices

[~, J] = size(income);

pre_slutsky = zeros(J, J, R);
for r=1:R
    numerator = exp(delta_blp_s9w10' + sigmaI_blp.*income_draws(r).*price_s9_w10 + sigmaB_blp.*v_draws(r).*branded_s9w10');
    share_r = numerator./(1+sum(numerator));
    pre_slutsky(:, :, r) = (alpha_blp + sigmaI_blp.* income_draws(r)).*price_s9_w10'.*(diag(share_r) - share_r*share_r');
end

slutsky_blp = (1./mkt_share_s9_w10).*mean(pre_slutsky,3);

%% Supply side


% we back out marginal costs by solving a linear system coming from firm's
% FOCs - formula is in the notes

% first we define an ownership matrix - omega_jk = 1 if products j and k
% belong to the same firm

Omega = [1 1 1 0 0 0 0 0 0 0 0;
         1 1 1 0 0 0 0 0 0 0 0;
         1 1 1 0 0 0 0 0 0 0 0;
         0 0 0 1 1 1 0 0 0 0 0;
         0 0 0 1 1 1 0 0 0 0 0;
         0 0 0 1 1 1 0 0 0 0 0;
         0 0 0 0 0 0 1 1 1 0 0;
         0 0 0 0 0 0 1 1 1 0 0;
         0 0 0 0 0 0 1 1 1 0 0;
         0 0 0 0 0 0 0 0 0 1 1;
         0 0 0 0 0 0 0 0 0 1 1];

wholesale_cost_s9_w10 = data_long.cost(data_long.store == 9 & data_long.week == 10);

% logit case
mc_logit = price_s9_w10 + inv(Omega.*slutsky_logit.*mkt_share_s9_w10./price_s9_w10)*mkt_share_s9_w10;

% blp case

mc_blp = price_s9_w10 + inv(Omega.*slutsky_blp.*mkt_share_s9_w10./price_s9_w10)*mkt_share_s9_w10;

%% MERGER

% New ownership matrix


Omega_merge = [1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         1 1 1 1 1 1 1 1 1 0 0;
         0 0 0 0 0 0 0 0 0 1 1;
         0 0 0 0 0 0 0 0 0 1 1];

% use estimated marginal costs and formula for logit elasticities to back out prices
% that solve firm's system of FOCs - use  previous prices as initial
% guesses

new_prices = fsolve(@(params)merger_logit(params(1:11), mc_logit, Omega_merge, delta_logit_s9w10, alpha_logit, price_s9_w10), price_s9_w10);
% test that with the status quo nothing changes - all equal!
test_prices= fsolve(@(params)merger_logit(params(1:11), mc_logit, Omega, delta_logit_s9w10, alpha_logit, price_s9_w10), price_s9_w10);







