function [omega, sj] = omegafinder(delta, sigma, price, nu, alpha)

ns = size(nu,2);
J = size(delta,1);

mu = sigma*price*nu;
expdelta_mu = exp(repmat(delta+ price*alpha,1,ns) + mu);
denominator = 1 + sum(expdelta_mu,1); 
indiv_shares = expdelta_mu./(repmat(denominator,J,1));
sj = 1./ns*sum(indiv_shares,2);


alpha = alpha + repmat(sigma*nu,J,1);

ind_share_alpha = sqrt(alpha).*indiv_shares;
sub_offdiag = ind_share_alpha*ind_share_alpha.';

ind_share_alpha = alpha.*indiv_shares;
sub_diag = diag(sum(ind_share_alpha,2));

omega = 1/ns*(sub_diag-sub_offdiag);

end