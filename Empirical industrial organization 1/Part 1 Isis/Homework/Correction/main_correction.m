clear all

data = csvread('data.csv');

year = data(:,1);

kt = size(year,1);

num_market = year - 2002;

T = max(num_market);

prodsMarket = accumarray(num_market,ones(kt,1));


marketStarts = zeros(T,1);
marketEnds = zeros(T,1);
marketStarts(1) = 1;
marketEnds(1) = prodsMarket(1);
for t = 2:T
    marketStarts(t) = marketEnds(t-1) + 1;
    marketEnds(t) = marketStarts(t) + prodsMarket(t) - 1;
end

num_group = data(:,2);

num_seg = data(:,3);

price = data(:,4);

xo = data(:,5:8);
%in order: cylinder, weight, HP, fuel cost

fuel_cost = data(:,8);

no_hh = data(:,9);

sj = data(:,10);

s0 = data(:,11);

iv = data(:,12:23);

sjg = zeros(kt,1);

for tt = 1:T
    Jt = marketStarts(tt): marketEnds(tt);
    sum_seg = accumarray(num_seg(Jt), sj(Jt));
    sjg(Jt) = sj(Jt)./sum_seg(num_seg(Jt));
end

l_sbar = data(:,24);

l_sbar_check = log(sjg);

max(abs(l_sbar-l_sbar_check))

y = data(:,25);

ye_fe = data(:,26:30);

brd_fe = data(:, 31:69);

xo_big = [ones(kt,1), xo, ye_fe, brd_fe];

z = [xo_big, iv];

n_z = size(z,2);

invZ = (z'*z)\eye(n_z);

%% Nested logit
% [theta_est, xi] = ivregress(y, xo_big, [price, l_sbar_check], z, invZ);
% 
% alpha = theta_est(end-1);
% 
% sigma = theta_est(end);
% 
% delta_hp = xo_big*theta_est(1:end-2) + xi;
% 
% delta_hp_check = y - alpha*price - sigma*l_sbar_check;
% 
% max(abs(delta_hp_check-delta_hp))
% 
% 
% sj_check = zeros(kt,1);
% for tt = 1:T
%     Jt = marketStarts(tt):marketEnds(tt);
%     sj_check(Jt) = shares_nested(delta_hp(Jt), price(Jt), alpha, sigma, num_seg(Jt));
% end
% 
% max(abs((sj_check-sj)./sj))
% 
% cost = zeros(kt,1);
% 
% for tt = 1:T
%     Jt = marketStarts(tt):marketEnds(tt);
%     kkt = size(Jt,2);
%     omeg = alpha * omega_nested(sigma, sj(Jt), sjg(Jt), num_seg(Jt));
%     owner = (repmat(num_group(Jt),1,kkt) == repmat(num_group(Jt)',kkt,1));
%     omeg_tile = omeg.*owner;
%     cost(Jt) = price(Jt) + omeg_tile\sj(Jt);
% end
% 
% min(cost)
% sum(cost<0)
% mean((price-cost)./price)
% mean(price./cost)
% 
% 
% opt = optimset('Display', 'iter',...
%     'TolFun',1e-4,'TolX',1e-6,...
%     'MaxFunEvals', 1000000,'MaxIter',5000,...
%     'LargeScale','off');
% 
% check = zeros(kt,1);
% for tt = 1:1
%      Jt = marketStarts(tt):marketEnds(tt);
%      price_check = fsolve(@(price_temp) foc(price_temp, cost(Jt), delta_hp(Jt), ...
%      num_group(Jt), alpha, sigma, num_seg(Jt)), price(Jt), opt);
%      check(Jt) = foc(price(Jt), cost(Jt), delta_hp(Jt), num_group(Jt), alpha, sigma, num_seg(Jt));
% end
% 
%%
rng(2020);
ns = 100;
nu_temp = randn(1, ns/2);
nu = [nu_temp,-nu_temp];

% nu = nu(abs(nu) < 2);
% N_temp = size(nu,2);
% ii = 0;
% while N_temp<ns
%     nu_temp = randn(1,ns - N_temp);
%     nu_new = [nu, nu_temp(abs(nu_temp)<2)];
%     N_temp = size(nu_new,2);
%     nu = nu_new;
%     ii = ii + 1
% end

grid = [-4:0.5:4];
nb_grid = size(grid,2);
fval_grid = zeros(nb_grid,1);
for ii = 1:nb_grid
   fval_grid(ii) = gmm_random(grid(ii), xo_big, price, price, nu,...
    sj, z, invZ, y, marketStarts, marketEnds);

end

figure
plot(grid, fval_grid)


[sigma_hat, fval] = fminunc(@(theta)gmm_random(theta, xo_big, price, price, nu,...
    sj, z, invZ, y, marketStarts, marketEnds), 1.5, opt);

mu_est = sigma_hat*price*nu;

delta_est = zeros(kt,1);
for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    delta_est(Jt) = inv_sh_by_market(sj(Jt), y(Jt), mu_est(Jt,:));
end


[beta_hat, xi_hat] = ivregress(delta_est, xo_big, price, z, invZ);

alpha_hat = beta_hat(end);

beta_fc = beta_hat(5);


delta_hp = delta_est - alpha_hat*price;

alpha_i = alpha_hat + nu*sigma_hat;

figure
hist(alpha_i)
mean(alpha_i > 0)


for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    kkt = size(Jt,2);
    omeg = omega_rc(delta_est(Jt), price(Jt), num_group(Jt), alpha_hat, sigma_hat, nu);
    cost(Jt) = price(Jt) + omeg\sj(Jt);
end

min(cost)
sum(cost<0)
mean((price-cost)./price)
mean(price./cost)

check = zeros(kt,1);
price_check = zeros(kt,1);
for tt = 1:T
     Jt = marketStarts(tt):marketEnds(tt);
     price_check(Jt) = fsolve(@(price_temp) foc_rc(price_temp, cost(Jt), delta_hp(Jt), ...
     num_group(Jt), alpha_hat, sigma_hat, nu), price(Jt), opt);
     check(Jt) = foc_rc(price(Jt), cost(Jt), delta_hp(Jt), num_group(Jt), alpha_hat, sigma_hat, nu);
end

max(abs(check))

max(abs(price_check-price))


Jt = marketStarts(T):marketEnds(T);

%% projet 1: merger between VW (21) & BMW (1)
num_group_bis = num_group(Jt);

num_group_bis(num_group(Jt) == 21) = 1;

sum(num_group_bis==21)

price_new1 = fsolve(@(price_temp) foc_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    num_group_bis, alpha_hat, sigma_hat, nu), price(Jt), opt);

max((price_new1-price(Jt))./price(Jt))
mean((price_new1-price(Jt))./price(Jt))
min((price_new1-price(Jt))./price(Jt))


grid = [0.01:0.005:0.05];
nb_grid = size(grid, 2);
csurplus = zeros(nb_grid,1);

for ii = 1:nb_grid

cost_new = cost(Jt);

cost_new(num_group_bis == 1) = (1-grid(ii)) * cost_new(num_group_bis == 1);

price_new11 = fsolve(@(price_temp) foc_rc(price_temp, cost_new, delta_hp(Jt), ...
    num_group_bis, alpha_hat, sigma_hat, nu), price(Jt), opt);

max((price_new11-price(Jt))./price(Jt))
mean((price_new11-price(Jt))./price(Jt))
min((price_new11-price(Jt))./price(Jt))

delta_new11 = delta_hp(Jt)+alpha_hat*price_new11;

csurplus(ii) = surplus_by_market(delta_new11, price_new11, sigma_hat, nu, alpha_hat);

end

csurplus_init = surplus_by_market(delta_est(Jt), price(Jt), sigma_hat, nu, alpha_hat);

figure
plot(grid, csurplus)
hold on
plot([min(grid);max(grid)], [csurplus_init;csurplus_init])
axis tight
hold off


% %% projet 2: merger between RENAULT (15) & PSA(14)
% num_group_bis = num_group(Jt);
% 
% num_group_bis(num_group(Jt) == 15) = 14;
% 
% sum(num_group_bis==15)
% 
% price_new2 = fsolve(@(price_temp) foc(price_temp, cost(Jt), delta_hp(Jt), ...
%     num_group_bis, alpha, sigma, num_seg(Jt)), price(Jt), opt);
% 
% max((price_new2 - price(Jt))./price(Jt))
% mean((price_new2 - price(Jt))./price(Jt))
% min((price_new2-price(Jt))./price(Jt))

%max price increase = 9.68%

% %% projet 3 + 4: effect of BMW & VW REN & PSA merging
% num_group_bis = num_group(Jt);
% 
% num_group_bis(num_group(Jt) == 15) = 14;
% 
% num_group_bis(num_group(Jt) == 21) = 1;
% 
% price_new3 = fsolve(@(price_temp) foc(price_temp, cost(Jt), delta_hp(Jt), ...
%     num_group_bis, alpha, sigma, num_seg(Jt)), price(Jt), opt);
% 
% max((price_new3-price(Jt))./price(Jt))
% mean((price_new3-price(Jt))./price(Jt))
% min((price_new3-price(Jt))./price(Jt))

% max price increase = 9.71
%% projet 5 + 6: effect of cross participations 1) REN (15) buys 30% PSA (14)
%i.e. REN internalizes 30% of PSA's profits
kkt = size(Jt,2);
kappa_mat = zeros(kkt,kkt);

kappa_mat(num_group(Jt)==15, num_group(Jt)==14) = 0.3;

mat_own = repmat(num_group(Jt),1,kkt)==repmat(num_group(Jt)',kkt,1);

owner = mat_own + kappa_mat;

price_new4 = fsolve(@(price_temp) foc_kappa_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    alpha_hat, sigma_hat, nu, owner), price(Jt), opt);

max((price_new4-price(Jt))./price(Jt))
mean((price_new4-price(Jt))./price(Jt))
min((price_new4-price(Jt))./price(Jt))


% PSA buys 30% REN in exchange

kappa_mat(num_group(Jt)==14, num_group(Jt)==15) = 0.3;

owner22 = mat_own + kappa_mat;

price_new5 = fsolve(@(price_temp) foc_kappa_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    alpha_hat, sigma_hat, nu, owner22), price_new4, opt);

max((price_new5 - price(Jt))./price(Jt))
mean((price_new5 - price(Jt))./price(Jt))
min((price_new5 - price(Jt))./price(Jt))

%without buy of Renault
kappa_mat = zeros(kkt,kkt);
kappa_mat(num_group(Jt)==14, num_group(Jt)==15) = 0.3;

owner21 = mat_own + kappa_mat;

price_new5b = fsolve(@(price_temp) foc_kappa_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    alpha_hat, sigma_hat, nu, owner21), price_new4, opt);

max((price_new5b - price(Jt))./price(Jt))
mean((price_new5b - price(Jt))./price(Jt))
min((price_new5b - price(Jt))./price(Jt))



%% projet 6: effect of cross participations between  VW (21) & BMW (1)
%1)VW buys BMW
%2)BMW buys VW
% 
% kkt = size(Jt,2);
% kappa_mat = zeros(kkt,kkt);
% 
% kappa_mat(num_group(Jt)==21, num_group(Jt)==1) = 0.3;
% 
% mat_own = repmat(num_group(Jt),1,kkt)==repmat(num_group(Jt)',kkt,1);
% 
% owner = mat_own + kappa_mat;
% 
% price_new6 = fsolve(@(price_temp) foc_kappa(price_temp, cost(Jt), delta_hp(Jt), ...
%     owner, alpha, sigma, num_seg(Jt)), price(Jt), opt);
% 
% max((price_new6-price(Jt))./price(Jt))
% mean((price_new6-price(Jt))./price(Jt))
% min((price_new6-price(Jt))./price(Jt))
% 
% 
% % PSA buys 30% REN in exchange
% 
% kappa_mat(num_group(Jt)==1, num_group(Jt)==21) = 0.3;
% 
% owner = mat_own + kappa_mat;
% 
% price_new7 = fsolve(@(price_temp) foc_kappa(price_temp, cost(Jt), delta_hp(Jt), ...
%     owner, alpha, sigma, num_seg(Jt)), price(Jt), opt);
% 
% max((price_new7-price(Jt))./price(Jt))
% mean((price_new7-price(Jt))./price(Jt))
% min((price_new7-price(Jt))./price(Jt))


