function [omega, sj] = omega_share(delta_pure,alpha,price,ns)

numerator = exp(delta_pure - alpha.*price); 
denominator = 1+sum(numerator,1); 
s_ij = numerator./denominator;
sj = (1/ns)*sum(s_ij,2);

right = - alpha.*s_ij*s_ij';
left = diag(sum(- alpha.*s_ij,2)); 
omega = (1/ns)*(left-right);






