function cs = consumer2(delta_np,price, beta, tax, alpha, fuel_cost,sigma,nu, ns )

a = exp(delta_np + beta *(price + tax)  + sigma*fuel_cost*nu);
b = sum(a,1); %sum of products for each person i
c = log(1+b);
alphai = alpha + sigma*nu; % alpha for each i
d = c./abs(alphai); %CS for each i
cs = sum(d)/ns;
