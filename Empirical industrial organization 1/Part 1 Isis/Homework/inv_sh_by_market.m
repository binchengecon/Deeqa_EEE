function res = inv_sh_by_market(sj, delta_init, mu)

kt = size(delta_init,1);

tol = 1e-7;

iter = 0;
deltaold = zeros(kt,1);
delta = delta_init;

while(max(abs(delta - deltaold)) > tol && iter < 5000)
    iter = iter+1;
    deltaold = delta;
    delta = delta + log(sj) - log(agg_sh_by_market(delta, mu));
end

% disp(['no iteration for inversion: ', num2str(iter)])

res = delta;