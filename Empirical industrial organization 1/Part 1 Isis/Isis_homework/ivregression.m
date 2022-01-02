function [beta, residual, var_matrix] = ivregression(y, xo, price, z, w)

beta_first = w * z' * price;

price_hat = z * beta_first;

x_hat = [xo, price_hat];

beta = (x_hat' * x_hat)\(x_hat' * y);

residual = y - [xo, price] * beta;

P = z*inv(z'*z)*z';
var_residual = 1/size(xo,1)*sum(residual.^2);
x_and_price = [xo, price];
var_matrix = var_residual*inv(x_and_price'*P*x_and_price);

