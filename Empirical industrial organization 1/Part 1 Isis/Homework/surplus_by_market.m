function [cs] = surplus_by_market(delta, price, sigma, nu, alpha)

kt = size(delta,1);
ns = size(nu,2);

mu =  price * sigma * nu;

numer = exp(repmat(delta,1,ns) + mu);

sum1 = 1 + sum(numer);

alpha_i = -(alpha + sigma * nu);

cs_i = log(sum1) ./  alpha_i;

cs = sum(cs_i*10000)/ns;