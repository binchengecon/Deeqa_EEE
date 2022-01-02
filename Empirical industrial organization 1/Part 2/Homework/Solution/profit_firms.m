function [s_hat_j, profit] = profit_firms(delta_np,price,alpha,I,mc)

numer = exp(delta_np - alpha'.*price); % 4 x 10000

denom = 1+sum(numer,1); % 1 x 10000

s_hat_ij = numer./denom; %  4 x 10000

s_hat_j= mean(s_hat_ij,2);  

profit = (price-mc).*s_hat_j*I;
end
