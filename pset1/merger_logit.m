function foc = merger_logit(new_price, mc, Omega, delta_old, alpha_logit, old_price)
%little trick just so we can express stuff in terms of old deltas, not on
%terms of observables X
numerator = exp(delta_old + alpha_logit*new_price - alpha_logit*old_price);
share= numerator./(1+sum(numerator));
new_slutsky = (1./share).*alpha_logit.*new_price'.*(diag(share) - share*share');
foc = share + (Omega.*new_slutsky.*share./new_price)*(new_price - mc);


end

