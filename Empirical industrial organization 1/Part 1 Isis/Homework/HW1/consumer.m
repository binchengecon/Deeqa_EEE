function cs = consumer(delta_np,alpha,price,sigma,nu,ns)

a = exp(delta_np + alpha*price + sigma*nu.*price);
b = sum(a,1); %sum of products for each person i
c = log(1+b);
alphai = alpha + sigma*nu; % alpha for each i
d = c./abs(alphai); %CS for each i
cs = sum(d)/ns;
