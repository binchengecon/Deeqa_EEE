function [delta_new, pb_cv] = blp_contraction(sigma, sj, price, delta_init, nu)
%is valid for inversion market by market

J = size(delta_init,1);
delta_old = zeros(J,1);
delta_new = delta_init;
iter = 0;

while max(abs(delta_new - delta_old))> 1e-7 && iter < 1000
    iter = iter + 1;
    delta_old = delta_new;
    delta_new = delta_old + log(sj) - log(shares(delta_old, sigma, price, nu));
end

pb_cv = 0;
if iter == 1000
    pb_cv = 1;
end


% disp(['no iter =', num2str(iter)])