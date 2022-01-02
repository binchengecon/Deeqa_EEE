function  res = gmm_blp(sigma, sj, x0 , x1 , endo, z, w, nu, delta_init, yearStarts, yearEnds)

% x0: exogenous variable
% x1: variable with random coefficient
% endo: endogenous variable
%       size [kt,2]
kt = size(sj,1);
T = size(yearStarts,1);
delta = zeros(kt,1);
pb_cv = zeros(T,1);

for tt = 1:T
Jt = yearStarts(tt):yearEnds(tt);    
[delta(Jt), pb_cv(tt)] = blp_contraction(sigma, sj(Jt), delta_init(Jt), nu, x1(Jt));

% endo size [kt,2]
end


if sum(pb_cv) == 0
[~, xi] = ivregression(delta, x0, endo, z);
moments = z' * xi; %dim nb_z *1
res = moments'*w*moments;
else 
    res  = 1e5;
end