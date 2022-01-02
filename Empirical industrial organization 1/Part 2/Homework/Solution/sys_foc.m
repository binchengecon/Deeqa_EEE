function res = sys_foc(price,delta_np,alpha,ns,own,mc)


numer = exp(delta_np - alpha'.*price); % 4 x 10000

denom = 1+sum(numer,1); % 1 x 10000

s_hat_ij = numer./denom; %  4 x 10000

s_hat_j= mean(s_hat_ij,2);  

Temp1 = -diag(sum(alpha'.*s_hat_ij,2)); 
Temp2 = alpha'.*s_hat_ij*s_hat_ij';

omega = (1/ns)*(Temp1+Temp2);  % J X J matrix


res = 1e10*(s_hat_j + omega.*own * (price-mc)); %objective function