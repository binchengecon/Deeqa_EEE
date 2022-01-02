function [beta, residual] = ivregression(y, xo, price, z)



beta_first = inv(z'*z) * z' * price;
price_hat = z * beta_first;
x_hat = [xo, price_hat];
beta = inv(x_hat' * x_hat) * (x_hat' * y);
residual = y - [xo, price] * beta;