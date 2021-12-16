function res = gmm_obj2(theta, price, qty, xd, xs, z, w, nfirm)

n_z = size(z, 2);

n_xd = size(xd,2);

beta_d = theta(1:n_xd);

beta_s = theta(n_xd+1:end);

b0 = beta_d(end);

resid_d = price - xd*beta_d;

resid_s = price - b0*qty./nfirm - xs*beta_s;

moment_s = sum(repmat(resid_s, 1, n_z).*z);

moment_d = sum(repmat(resid_d, 1, n_z).*z);

moment = [moment_d, moment_s];

res = moment*w*moment';