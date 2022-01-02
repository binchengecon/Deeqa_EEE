function residual = help_RPM_foc(delta_pure, alpha, price, ns, mc)

[omega, sj]= omega_share(delta_pure,alpha,price,ns);
residual = 1e5*(sj + omega'* (price-mc));

