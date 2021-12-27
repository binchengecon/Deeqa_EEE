%this function calculates the estimator for a linear equation using
%instrumental variablea
%It outputs the estimator and the residuals
function [beta_second, residuals] = ivregress(y, xe, xo, z)
%the argument of the function are the dependent variable y, the endogenous
%variable(s) xe, the exogenous explanatory variables xo and the
%excluded instruments z

%construct the full matrix of instruments from xo and z
zz = [xo, z];

%calculate the first stage estimator: (Z'Z)^(-1)*Z'*xe
beta_first = (zz'*zz)\(zz'*xe);

%calculate xe_hat = Z*beta_first
xe_hat = zz * beta_first ;

%construct the matrix of regressors for the second stage
x_hat = [xo, xe_hat];

%calculate the second stage estimator
beta_second = (x_hat' * x_hat) \( x_hat' * y);

%calculate the residuals
residuals = y - [xo, xe]*beta_second; 