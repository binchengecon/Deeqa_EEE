function omeg = omega(delta, price, num_group, alpha, sigma, nu)

kt = size(num_group,1);
ns = size(nu,2);

% i represent draws
% j represent product

alpha_i = alpha + sigma*nu;

mu = sigma*price*nu;

numer = exp(repmat(delta,1,ns) + mu); 

denom = repmat(1 + sum(numer),kt,1); 

sij = numer./(denom);   

% ownership matrix
owner = (repmat(num_group,1,kt) == repmat(num_group',kt,1));

% duplicate vector by group for multiple loop

alpha_sij = repmat(alpha_i,kt,1) .* sij;

sum_alpha_sij = sum(alpha_sij,2);

omeg_tilde = (1/ns) * (diag(sum_alpha_sij) - (alpha_sij * sij'));

omeg = omeg_tilde .* owner; %omega negative on diag