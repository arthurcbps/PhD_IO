% GMM objective function

% we will use out the Nevo 'trick' of using the focs to express the linear
% parameters (alpha, beta) as a function of the linear ones, which allows
% us to minimize the GMM criterion function considering only the non linear
% ones

% X1 as product characteristics - brand dummies, promotion, price
% Z are instruments to construct the weighting matrix
function moment = gmm(sigmaI, sigmaB, X1, Z, mkt_share, prices, income, branded, R)


delta_wide = iterate_delta(sigmaI, sigmaB, mkt_share, prices, income, branded, R);

% use Nevo's tip of setting the function to be arbitrarily high if mean
% utility is ill-defined
if max(isnan(delta_wide))== 1
    moment = 1e10;
else
    delta_long = delta_wide(:);
    %optimal weighting matrix
    W = inv(Z'*Z);
    % expression for linear parameters (the trick)

    lin_par = inv(X1'*Z*W*Z'*X1)*X1'*Z*W*Z'*delta_long;

    % back out xi's to form moment conditions

    xi = delta_long - X1*lin_par;

    moment = xi'*Z*W*Z'*xi;

end

