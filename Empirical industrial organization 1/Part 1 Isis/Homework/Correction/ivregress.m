function [beta, residual, std] = ivregress(y, xo, price, z, invZ)

[kt, n_endo]  = size(price);
K = size(xo,2) + size(price,2);


price_hat = zeros(kt,n_endo);
for ii=1:n_endo
beta_first = invZ * z' * price(:,ii);

price_hat(:,ii) = z * beta_first;
end 
x_hat = [xo, price_hat];

beta = (x_hat' * x_hat) \ (x_hat' * y);

residual = y - [xo, price] * beta;

RSS = sum(residual.^2);

m = ((repmat(residual.^2,1,K) .* x_hat)' * x_hat);

varcov = inv(x_hat'*x_hat) * m * inv(x_hat'*x_hat);

% int = ([x_hat]' * [x_hat]);
% 
% varcov = s2 * inv(int);


std = sqrt(diag(varcov));