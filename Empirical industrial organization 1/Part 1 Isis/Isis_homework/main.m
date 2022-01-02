%%Version of Hoang Phan

clear all
clc
%%
cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Isis_homework'
data = csvread("base_project_deeqa_notext.csv");

sj = data(:,10);

xo_temp = data(:,5:8);

iv = data (:,11:22);

price = data(:,4);

num_firm = data(:,2);

num_brand = data(:,3);

num_market = data(:,1) - 2002; %transform year to market index

kt = size(sj,1);

xo = [ones(kt,1), xo_temp]; %add intercept to the xo

z = [xo, iv];

w = inv(z'*z); %avoid computing it during the estimation

%% Dummy variables
brand_dummies = dummyvar(num_brand);
time_dummies = dummyvar(num_market);
%remove the first column to avoid collinearity
brand_dummies(:,1)=[];
time_dummies(:,1) =[];
dummies = [brand_dummies, time_dummies];
%% Useful index: locate begin and end of each market

T = max(num_market);
prodsMarket = zeros(T,1);
kt = size(num_market,1);
for t = 1:T
    prodsMarket(t) = sum(num_market==t);
end

marketStarts = zeros(T,1);
marketEnds = zeros(T,1);
marketStarts(1) = 1;
marketEnds(1) = prodsMarket(1);

for t=2:T
    marketStarts(t) = marketEnds(t-1) + 1;
    marketEnds(t) = marketStarts(t) + prodsMarket(t) - 1;
end

%outside option
s0 = zeros(kt,1);
for t=1:T
    s0_temp = 1 - sum(sj(marketStarts(t):marketEnds(t)));
    s0(marketStarts(t):marketEnds(t)) = s0_temp;
end

y = log(sj./s0);


%% Inversion of the market share
%first, create the function that computes market share for delta, mu
%draw 100 normal random draws to approximate the integral
ns = 100;

%draw in the normal(0.1)
rng(2020);
nu = randn(1, ns);

delta_init = y;

delta_new = zeros(kt,1);
for tt = 1: T
    Jt = marketStarts(tt):marketEnds(tt);
    delta_new(Jt) = blp_contraction(0.4, sj(Jt), price(Jt), delta_init(Jt), nu);
end

%% Question 1.1
% Estimation of the RC model

opt = optimset('display', 'iter', 'tolX', 1e-4);

% [gmm_init,~] = gmm_blp(0.1, sj, xo, price, dummies, z, invZ, nu, delta_init, ...
%                         marketStarts, marketEnds, false, true);
% 
% gmm_fun = zeros(11,1);
% jj = 0;
% for ii = 0:0.1:1
%     jj = jj +1;
%     [gmm_fun(jj),~] = gmm_blp(ii, sj, xo, price, dummies, z, invZ, nu, delta_init, ...
%                                 marketStarts, marketEnds, false, true);
% end
% figure
% plot(0:0.1:1, gmm_fun)

theta_hat = fminunc(@(theta)gmm_blp(theta, sj, xo, price,z, w, nu, delta_init, ...
    marketStarts, marketEnds), 0.1, opt);

[gmm_min,beta,var_matrix] = gmm_output(theta_hat, sj, xo, price, z, w, nu, delta_init,...
                        marketStarts, marketEnds);
% Coefficient: constant, cylinder, weight, horsepower, fuel cost, price
beta
% Asymptotic variance:
diag(var_matrix)
% P-value:
T_stat = beta./sqrt(diag(var_matrix));
p_value = 1-2.*abs(0.5 - normcdf(T_stat))

%% Estimation of the Fixed Effects model

opt = optimset('display', 'iter', 'tolX', 1e-4);
%2SLS with dummies
xo_fe = [xo, dummies];
z_fe = [xo, dummies, iv];
w_fe = inv(z_fe'*z_fe);
% gmm_init = gmm_blp(0.1, sj, xo_fe, price, dummies, z_fe, w_fe, nu, delta_init, ...
%                         marketStarts, marketEnds, true, true);
% sigma_range = 0:0.1:4;
% n_sigma = length(sigma_range);
% gmm_fun = zeros(n_sigma,1);
% jj = 0;
% for ii = sigma_range
%     jj = jj +1;
%     gmm_fun(jj) = gmm_blp(ii, sj, xo_fe, price, dummies, z_fe, w_fe, nu, delta_init, ...
%                                 marketStarts, marketEnds, true, true);
% end
% figure
% plot(sigma_range, gmm_fun)

theta_hat_feiv = fminunc(@(theta)gmm_blp(theta, sj, xo_fe, price, z_fe, w_fe, nu, delta_init, ...
    marketStarts, marketEnds), 0.1, opt);

[gmm_min_feiv,beta_feiv,var_matrix_feiv, delta] = gmm_output(theta_hat_feiv, sj, xo_fe, price, z_fe, w_fe, nu, delta_init,...
                        marketStarts, marketEnds);

beta_feiv_nodum = beta_feiv([1:5,size(beta_feiv,1)])
% Asymptotic variance:
vardiag = diag(var_matrix_feiv);
coef_var = vardiag([1:5,size(beta_feiv,1)]);
% P-value:
T_stat = beta_feiv_nodum./sqrt(coef_var);
p_value = 1-2.*abs(0.5 - normcdf(T_stat))

%% Estimation with differentiation IV
% Isolation local differentiation IV

beta_first = w_fe * z_fe' * price;

price_hat = z_fe * beta_first;

char = [xo, price_hat]; %characteristics including price_hat
char(:,1) = []; %remove constant variable
sd_vector = std(char);
diffIV = zeros(size(char));

for t = 1:T
    char_t = char(marketStarts(t):marketEnds(t),:);
    no_products = size(char_t,1);
    no_char = size(char_t,2);
    
    for k = 1:no_char
        tempMatrix = repmat(char_t(:,k),1,no_products);
        diff = tempMatrix' - tempMatrix; 
        absdiff = abs(diff);
        indicator = absdiff < sd_vector(k);
        diffIV_k = sum(indicator.*absdiff);
        diffIV(marketStarts(t):marketEnds(t),k) = diffIV_k';  
    end
end
    
%Estimation with diffIV

opt = optimset('display', 'iter', 'tolX', 1e-4);

newIV = [z_fe, diffIV];
inv_newIV = inv(newIV'*newIV);
gmm_init = gmm_blp(0.1, sj, xo_fe, price, newIV, inv_newIV, nu, delta_init, ...
                        marketStarts, marketEnds);
sigma_range = 0:0.1:4;
n_sigma = length(sigma_range);
gmm_fun = zeros(n_sigma,1);
jj = 0;
for ii = sigma_range
    jj = jj +1;
    gmm_fun(jj) = gmm_blp(ii, sj, xo_fe, price, newIV, inv_newIV, nu, delta_init, ...
                                marketStarts, marketEnds);
end
figure
plot(sigma_range, gmm_fun)

theta_hat_diffiv = fminunc(@(theta)gmm_blp(theta, sj, xo_fe, price, newIV, inv_newIV, nu, delta_init, ...
                                marketStarts, marketEnds), 0.1, opt);

[gmm_min_diffiv,beta_diffiv, var_matrix_diffiv] = gmm_output(theta_hat_diffiv, sj, xo_fe, price, newIV, inv_newIV, nu, delta_init, ...
                                marketStarts, marketEnds);
    
beta_diffiv_nodum = beta_diffiv([1:5,size(beta_diffiv,1)])
% Asymptotic variance:
vardiag = diag(var_matrix_diffiv);
coef_var = vardiag([1:5,size(beta_diffiv,1)]);
% P-value:
T_stat = beta_diffiv_nodum./sqrt(coef_var);
p_value = 1-2.*abs(0.5 - normcdf(T_stat))
    