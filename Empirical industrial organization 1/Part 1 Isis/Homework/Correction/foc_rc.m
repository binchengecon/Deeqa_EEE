function res = foc_rc(price_temp, cost, delta_hp, num_group, alpha, sigma, nu)

kkt = size(num_group,1);
ns = size(nu,2);

alpha_i = alpha + sigma*nu;

delta = delta_hp + alpha*price_temp;

mu = sigma*price_temp*nu;

[sj, sij] = agg_sh_by_market(delta, mu);

owner = (repmat(num_group,1,kkt) == repmat(num_group',kkt,1));

alpha_sij = repmat(alpha_i,kkt,1).* sij ;

sum_alpha_sij = sum(alpha_sij,2);

omeg_tilde = (1/ns) * (diag(sum_alpha_sij)-(alpha_sij * sij'));

omeg = omeg_tilde .* owner;

res = 1e6 * (omeg * (price_temp - cost) + sj) ;




