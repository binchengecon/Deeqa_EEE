function residual = help_foc(delta_pure, alpha, price, ns, own_matrix, mc)

[omega, sj]= omega_share(delta_pure,alpha,price,ns);
residual = 1e5*(sj + omega.*own_matrix * (price-mc));

