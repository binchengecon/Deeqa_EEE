function res = sys_foc2(price,delta_np,alpha,I,mc)

[omega, sj]= mega(delta_np,price,alpha,I);
res = 1e10*(sj + omega' * (price-mc)); %objective function

