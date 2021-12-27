function res = sys_foc(price,delta_np,alpha,I,own,mc)

[omega, sj]= mega(delta_np,price,alpha,I);
res = 1e10*(sj + omega.*own * (price-mc)); %objective function