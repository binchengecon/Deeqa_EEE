function res = sys_foc(price,delta,nu,sigma,alpha,ns,own,mc)

[omega, sj]= mega(delta,nu,price,sigma,alpha,ns);
res = 1e10*(sj + omega.*own * (price-mc)); %objective function