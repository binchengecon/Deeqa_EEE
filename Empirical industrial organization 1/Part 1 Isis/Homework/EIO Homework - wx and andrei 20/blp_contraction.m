function [delta_new, pb_conv] = blp_contraction(sigma, sj, price, delta_init, nu)

J         = size(delta_init,1);
delta_new = delta_init;
delta_old = zeros(J,1);
iter      = 1;
err       = 1e-7;

while max(abs(delta_new-delta_old)) > err && iter < 1000
    delta_old = delta_new;
    delta_new = delta_old + log(sj) - log(shares(delta_old, sigma, price, nu));
    iter      = iter +1;
end

pb_conv = 0;
if iter == 999
   pb_conv = 1; 
end

end