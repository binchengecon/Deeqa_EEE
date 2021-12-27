function aggregate_market_share = shares(delta, sigma, p, nu)
ns = size(nu,2);
J = size(delta,1);

%nu = vector column with draws, size (1, ns)

mu = sigma*p*nu;%matrix size J, ns

expdelta_mu = exp(repmat(delta,1, ns) + mu);

denominator = 1 + sum(expdelta_mu); %size 1, ns

indiv_shares = expdelta_mu./(repmat(denominator, J, 1)); %size J, ns

aggregate_market_share = 1./ns*sum(indiv_shares, 2); %vector J,1