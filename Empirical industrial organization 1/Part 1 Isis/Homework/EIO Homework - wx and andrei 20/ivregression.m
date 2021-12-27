function [beta, residual] = ivregression(y, x0, price, z)

% First stage regression
beta_first = inv(z'*z)*z'*price;
price_hat = z*beta_first;

% Second stage regression
x_hat = [x0, price_hat];
beta = inv(x_hat'*x_hat)* (x_hat'*y);
residual = y - [x0, price] * beta;

end