function aggregate_market_share = shares(delta, sigma, p, nu)

ns = size(nu,2);
J = size(delta,1);

mu = sigma*p*nu; % p in rows, nu in columns, columns consists of  

expdelta_mu = exp(repmat(delta,1,ns) + mu);
% Given a particular draw for an individual, calculate his total option
% value.
denominator = 1 + sum(expdelta_mu,1); 

indiv_shares = expdelta_mu./(repmat(denominator,J,1));

aggregate_market_share = 1./ns*sum(indiv_shares,2);

end