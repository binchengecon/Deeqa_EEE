function [res] = gmm_random(sigma, xo, x1, price, nu, sj, z, invZ, delta_init, marketStarts, marketEnds)

mu = sigma*x1*nu;

kt = size(sj,1);

T = size(marketStarts);

delta = zeros(kt,1);
for t = 1:T
    Jt = marketStarts(t):marketEnds(t);
    delta(Jt) = inv_sh_by_market(sj(Jt), delta_init(Jt), mu(Jt,:));
end


[beta, xi] = ivregress(delta, xo, price, z, invZ);

res = xi'*z*invZ*z'*xi;


