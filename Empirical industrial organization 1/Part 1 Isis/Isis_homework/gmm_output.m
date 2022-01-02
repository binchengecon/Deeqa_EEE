function [res,beta,var_matrix,delta] = gmm_output(sigma, sj, xo, price, z, w, nu,delta_init, marketStarts, marketEnds)

kt = size(sj,1);
T = size(marketStarts,1);

delta = zeros(kt,1);
pb_cv = zeros(T,1);

for tt = 1:T
    
Jt = marketStarts(tt):marketEnds(tt);    
[delta(Jt), pb_cv(tt)] = blp_contraction(sigma, sj(Jt), price(Jt), delta_init(Jt), nu);

end
if sum(pb_cv) == 0

    [beta, xi, var_matrix] = ivregression(delta, xo, price, z, w);
   
    moments = z' * xi; %dim nb_z *1
    res = moments'*w*moments;
else 
    res  = 1e5;
end