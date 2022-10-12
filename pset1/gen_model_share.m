function model_share = gen_model_share(delta, sigmaI, sigmaB, prices, income, branded)


% first step - fix non linear parameters sigmaB, sigmaI and draw values to
% simulate market shares
% we need to draw R times the individual-level terms in the nonlinear component of utility  mu_ij -
% those are the Gaussian error term and Income.
% we also take as given the linear term - delta

R = 50;
% normal shocks and income vary by individual only, so we draw R times for
% each

[nt, J] = size(prices);
% this ensures we draw the same numbers on every step of the contraction
% mapping
rng('default')

v = randn(nt, R);

% we draw income from the empirical distribution, assuming it is iid across
% observations
indices = randi(length(income), nt, R);

rand_income = income(indices);
% now we compute model-implied market shares for each combination of
% product j, week-store n, draw R


share_r = zeros(nt, J, R);
for r=1:R
    numerator = exp(delta + sigmaI.*rand_income(:, r).*prices + sigmaB.*v(:, r).*branded);
    share_r(:,:,r) = numerator./(1+sum(numerator, 2));
end

model_share = mean(share_r,3);

end