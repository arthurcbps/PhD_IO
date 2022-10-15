
%%% Now we used our implied model shares and the observed ones to iterate on a contraction and
%%% back out mean utility delta
function delta_new = iterate_delta(sigmaI, sigmaB, mkt_share, prices, income, branded, R)

% simulate market shares
% we need to draw R times the individual-level terms in the nonlinear component of utility  mu_ij -
% those are the Gaussian error term and Income.

% normal shocks and income vary by individual only, so we draw R times for
% each

nt = size(prices, 1);
v_draws = randn(nt, R);

% we draw income from the empirical distribution, assuming it is iid across
% observations
% set a seed for reproducibility

rng(1)
indices = randi(length(income), nt, R);

income_draws = income(indices);

delta_new = 0.1*ones(size(prices));

count = 1;
agg_error = 1;

while (agg_error > 1e-13) && count <= 10000
    delta_old = delta_new;
    delta_new = delta_old + log(mkt_share) - ...
        log(gen_model_share(delta_old, sigmaI, sigmaB, prices, branded, income_draws, v_draws, R));
    count= count+1;
    error = abs(delta_old - delta_new);
    agg_error = max(error, [], 'all');
disp(num2str(agg_error))
end

end



