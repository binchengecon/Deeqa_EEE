
function residual = CounterFactualSystem(price_temp, cost, co2, tax_rate, fuel_cost, delta_np, num_manuf, beta, sigma, nu)

kt = size(num_manuf,1);
ns = size(nu,2);
tax = co2*tax_rate/10000;
numer = exp(delta_np + beta*(price_temp + tax)  + sigma*fuel_cost*nu);
denom = 1+sum(numer,1);
% i for draws
s_hat_ij = numer./denom;

s_hat_j= mean(s_hat_ij,2);  

Temp_omega1 = diag(sum(beta.*s_hat_ij,2));
Temp_omega2 = beta.*s_hat_ij*s_hat_ij';

omega_tilde = (1/ns)*(Temp_omega1-Temp_omega2);

own = (repmat(num_manuf, 1, kt) == repmat(num_manuf', kt, 1));

omega = omega_tilde.*own;

residual = 1000000 * (omega * (price_temp - cost) + s_hat_j) ;
