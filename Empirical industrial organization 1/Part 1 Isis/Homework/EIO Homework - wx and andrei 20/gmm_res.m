function res = gmm_res(theta, sj, X, p, Z, invZ, nu, delta_init, marketStarts, marketEnds)

    [res,~,~] = gmm_blp(theta, sj, X, p, Z, invZ, nu, delta_init, marketStarts, marketEnds);

end