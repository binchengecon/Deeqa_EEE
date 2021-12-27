function [beta, residual] = regress(y, x)

beta = (x'*x) \ (x' * y);

residual = y - x * beta;