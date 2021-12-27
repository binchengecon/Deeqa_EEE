function [omega, sj] = mega(delta,nu,price,sigma,alpha,ns)

up = exp(delta + alpha*price + sigma*nu.*price);
down = 1+sum(up,1);
s_hat_ij = up./down;

second = (alpha+sigma*nu).*s_hat_ij*s_hat_ij';
first = diag(sum((alpha+sigma*nu).*s_hat_ij,2));



sj = (1/ns)*sum(s_hat_ij,2);
omega = (1/ns)*(first-second);




