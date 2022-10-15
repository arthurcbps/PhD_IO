function model_share = gen_model_share(delta, sigmaI, sigmaB, prices, branded, income_draws, v_draws, R)


% fix non linear parameters sigmaB, sigmaI, and income and
% gaussian shock draws
% we compute model-implied market shares for each combination of
% product j, week-store n, draw r


% we also take as given the linear term - delta
[nt, J] = size(prices);

share_r = zeros(nt, J, R);
for r=1:R
    numerator = exp(delta + sigmaI.*income_draws(:, r).*prices + sigmaB.*v_draws(:, r).*branded);
    share_r(:,:,r) = numerator./(1+sum(numerator, 2));
end

model_share = mean(share_r,3);

end