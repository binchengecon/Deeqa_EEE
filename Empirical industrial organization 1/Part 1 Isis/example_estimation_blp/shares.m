%this function calculate the market shares for a given delta for the model
%with random coefficient on price ONLY

function res = shares(delta, sigma, p, nu)

ns = size(nu,2);
%number of draws

J = size(delta,1);

%rk: nu must be a vector column with draws, it is of size (1, ns)
mu = sigma*p*nu;
%mu = matrix of size (J, ns)

expdelta_mu = exp(repmat(delta,1, ns) + mu);

denominator = 1 + sum(expdelta_mu); %size (1, ns)

indiv_shares = expdelta_mu./(repmat(denominator, J, 1)); %size (J, ns)

res = 1./ns*sum(indiv_shares, 2); %vector (J,1)