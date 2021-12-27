function [omega, sj] = mega(delta_np,price,alpha,I)


up = exp(delta_np + alpha'.*price); % 5 x 2000
down = 1+sum(up,1); % 1 x 2000, col sum
s_hat_ij = up./down; %  5 x 2000

second = alpha'.*s_hat_ij*s_hat_ij';
first = diag(sum(alpha'.*s_hat_ij,2)); % row sum

sj = (1/I)*sum(s_hat_ij,2);
omega = (1/I)*(first-second);  % J X J matrix




