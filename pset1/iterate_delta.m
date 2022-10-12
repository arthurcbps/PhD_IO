
%%% Now we used our implied model shares and the observed ones to iterate on a contraction and
%%% back out mean utility delta
function delta_new = iterate_delta(sigmaI, sigmaB, mkt_share, prices, income, branded)

delta_new = 0.1*ones(size(prices));

count = 1;
agg_error = 1;

while (agg_error > 0.001) && count <= 10000
    delta_old = delta_new;
    delta_new = delta_old + log(mkt_share) - log(gen_model_share(delta_old, sigmaI, sigmaB, prices, income, branded));
    count= count+1;
    error = abs(delta_old - delta_new);
    agg_error = max(error, [], 'all');
disp(num2str(agg_error))
end

end



