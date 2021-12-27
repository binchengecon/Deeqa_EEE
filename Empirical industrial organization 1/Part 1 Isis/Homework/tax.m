function tax_revenue = tax (price_temp, co2, tax_rate, fuel_cost, delta_np, beta, sigma, nu,nb_households)


tax = co2*tax_rate/10000;
numer = exp(delta_np + beta*(price_temp + tax)  + sigma*fuel_cost*nu);
denom = 1+sum(numer,1);
% i for draws
s_hat_ij = numer./denom;

s_hat_j= mean(s_hat_ij,2);  

q_hat_j = nb_households*s_hat_j;

tax_revenue = sum(tax.*q_hat_j)*10000/1000000;

