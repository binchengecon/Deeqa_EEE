function res = foc_kappa_rc(price_temp, cost, delta_hp, alpha, sigma, nu, owner)

kkt = size(delta_hp,1);
ns = size(nu,2);

alpha_i = alpha + sigma*nu;

delta = delta_hp + alpha*price_temp;

mu = sigma*price_temp*nu;

[sj, sij] = agg_sh_by_market(delta, mu);

alpha_sij = repmat(alpha_i,kkt,1).* sij ;

sum_alpha_sij = sum(alpha_sij,2);

omeg_tilde = (1/ns) * (diag(sum_alpha_sij)-(alpha_sij * sij'));

omeg = omeg_tilde .* owner;

res = 1e6 * (omeg * (price_temp - cost) + sj) ;




