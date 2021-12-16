%this function provides the systemn of first order conditions under a logit
%demand model
function res = sys_foc(price, delta_np, alpha, cost, ownership)

%compute mean product utilities including the price
edelta = exp(delta_np + alpha*price);

%compute denominator for market share equation
sum_edelta = 1+sum(edelta); %adds the outside good

%compute market shares
sj = edelta./sum_edelta;

%compute the matrix of partial derivatives ds_j/dp_k
omega = alpha * (diag(sj) - sj *sj');

%write de system of FOC we have to set to 0. Note that we multiply by a large value because
%otherwise the equation in the order of magnitude of sj---> very small!
res = 1e10*(sj + omega.*ownership * (price-cost));





