function res = sys_foc(price, delta, alpha, cost, ownership, nu, sigma)


[omega, sj] = omegafinder(delta, sigma, price, nu, alpha);

res = 1e5*(sj + omega.*ownership * (price-cost));
end