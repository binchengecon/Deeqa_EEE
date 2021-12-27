%% Change directory
cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Homework\EIO Homework - wx and andrei 20'
%% Read data
data  = readtable('base_project_deeqa.csv');

% Column names
names = data.Properties.VariableNames';

% Check the way market share is calculated using the outside option
% mktshare
num_obs = size(data,1);

data.s0 =  zeros(num_obs,1);

market_years = unique(data.year);
num_years = size(market_years,1);

% Create outside option as 26th column

for i = market_years'
    data_in_year = data(data.year==i,:);
    data(data.year==i,'s0') = table(repmat(1-sum(data_in_year.sj),sum(data.year==i),1));
end

% Indices for products in each year, for convenience

marketStarts = zeros(num_years,1);
marketEnds = zeros(num_years,1);
index = 1;
for i = market_years'
    marketStarts(index) = find(data.year==i,1,'first');
    marketEnds(index) = find(data.year==i,1,'last');
    index = index +1;
end

% Call out outside optiond
sj = data.sj;
s0 = data.s0;

y = log(sj./s0);


%% Q1: Random coefficient on price, linear in 5(6) variables, standard BLP instruments

characteristics = {'cylinder';'weight';'horsepower';'fuel_cost'};

market_time = unique(data.year); % years of markets

X = [table(ones(size(data,1),1)) data(:,characteristics)];
p = data(:,'price');
Z = [X data(:,14:25)];
 
p = p{:,:};
X = X{:,:};
Z = Z{:,:};

invZ = inv(Z'*Z);

rng(2020)
nu = randn(1,100); % Assume that the distribution of random coefficient is N(0,1), row of random numbers.

delta_init = y;
opt = optimset('display', 'iter', 'tolX', 1e-4);
[res_test, beta_test] = gmm_blp(0.5, sj, X, p, Z, invZ, nu, delta_init, marketStarts, marketEnds);

theta_hat = fminunc(@(theta)gmm_res(theta, sj, X, p, Z, invZ, nu, delta_init, ...
    marketStarts, marketEnds), 0.1, opt);

[res, beta] = gmm_blp(theta_hat, sj, X, p, Z, invZ, nu, delta_init, ...
    marketStarts, marketEnds);

gmm_fun = zeros(11,1);
jj = 0;
for ii = 0:0.1:1
    jj = jj +1;
    gmm_fun(jj) = gmm_blp(ii, sj, X, p, Z, invZ, nu, delta_init, marketStarts, marketEnds);
end
figure
plot(0:0.1:1, gmm_fun)

beta
theta_hat

%% Q2. Estimation of the RC with brand and time fixed effects.

% Generate time dummies

X = [X zeros(num_obs, num_years)];

year = min(data.year);
index = 6;
for i = 2003:max(data.year)
    X(:,index) = (data.year==year);
    year = year + 1;
    index = index + 1;
end
X(:,end) = [];
% Generate brand dummies

maxbrand = table2array(grpstats(data,{},@(x)size(unique(x),1),'DataVars',{'num_brand'}));
maxbrand = maxbrand(2);
X = [X zeros(num_obs, maxbrand)];
brand = 1;
index = index - 1;
for i = 1:maxbrand
    X(:,index) = (data.num_brand==brand);
    brand = brand + 1;
    index = index + 1;
end
X(:,end) = [];

Z = [X table2array(data(:,14:25))];
invZ = inv(Z'*Z);

delta_init = y;
opt = optimset('display', 'iter', 'tolX', 1e-4);
[res_test, beta_test] = gmm_blp(0.7, sj, X, p, Z, invZ, nu, delta_init, marketStarts, marketEnds);

gmm_fun = zeros(81,1);
jj = 0;
for ii = 0:0.05:4
    jj = jj +1;
    gmm_fun(jj) = gmm_blp(ii, sj, X, p, Z, invZ, nu, delta_init, marketStarts, marketEnds);
end
figure
plot(0:0.05:4, gmm_fun)

theta_hat = fminunc(@(theta)gmm_res(theta, sj, X, p, Z, invZ, nu, delta_init, ...
    marketStarts, marketEnds), 2.25, opt);

[res, beta, delta_cal] = gmm_blp(theta_hat, sj, X, p, Z, invZ, nu, delta_init, ...
    marketStarts, marketEnds);

beta2 = beta;
theta2 = theta_hat;

%% Q3. Differentiated IV

sd_instruments = zeros(num_years,4);
% 6 markets with 6 sets of 4 instrument sds
for i = 1:6
   for j = 1:4
       % location of instruments 8-cylinder : 11-fuel cost
       sd_instruments(i,j) = std(table2array(data(marketStarts(i):marketEnds(i),7+j)));
   end
end

diff_instrument = zeros(size(data,1),4);

for i = 1:6
    for j = 1:4
        current_market = table2array(data(marketStarts(i):marketEnds(i),7+j));
        for k = 1:size(current_market,1)
            diff_instrument(marketStarts(i)+k-1,j) = sum(abs(current_market(abs(current_market-current_market(k))<sd_instruments(i,j))-current_market(k)),1);
        end
    end
end

price_instrument_sd = zeros(num_years,1);
beta_price = inv(Z'*Z)*Z'*p;
price_hat = Z*beta_price;
price_instrument = zeros(size(data,1),1);
for i = 1:6
   price_instrument_sd(i) = std(price_hat(marketStarts(i):marketEnds(i)));  
end

for i = 1:6
   current_market = price_hat(marketStarts(i):marketEnds(i));
   for k = 1:size(current_market,1)
    price_instrument(marketStarts(i)+k-1) = sum(abs(current_market(abs(current_market-current_market(k))<price_instrument_sd(i))-current_market(k)),1);
   end
end

diffIV = [Z diff_instrument price_instrument];
invIV = inv(diffIV'*diffIV);
opt = optimset('display', 'iter', 'tolX', 1e-4);
theta_hat_iv = fminunc(@(theta)gmm_res(theta, sj, X, p, diffIV, invIV, nu, delta_init, ...
    marketStarts, marketEnds), 2.25, opt);

[res, beta_iv, delta_cal_iv] = gmm_blp(theta_hat_iv, sj, X, p, diffIV, invIV, nu, delta_init, ...
    marketStarts, marketEnds);

results = [beta2 beta_iv];

%% Counterfactual Simulations

market2008 = marketStarts(6):marketEnds(6);

shares2008 = sj(market2008);
prices2008 = p(market2008);

manu2008 = unique(data.num_manuf(data.year==2008));
ownershipmatrix1 = zeros(length(market2008));

X2008 = X(marketStarts(6):marketEnds(6),:);
manulist = data.num_manuf(market2008,:);
delta2008 = delta_cal(market2008);
alpha = beta2(end);
delta2008_np = delta2008 - prices2008*alpha;


[omega,sj_pred2008] = omegafinder(delta2008_np, theta2, prices2008, nu, alpha);
for i = 1:length(market2008)
    ownershipmatrix1(:,i) = (manulist == manulist(i));
end

o_tilde = ownershipmatrix1.*omega;
C2008 = prices2008+o_tilde\shares2008;
pcd2008 = prices2008-C2008;
num_hh2008 = unique(data.nb_households(data.year==2008));
profits2008 = pcd2008.*shares2008*num_hh2008;


a_pre_price = mean(prices2008);
a_pre_cost = mean(C2008);
a_pre_mp = mean(pcd2008./(C2008));
t_pre_ind = sum(profits2008,1);
a_pre_cs = consumer(delta2008_np,theta2,prices2008,nu,alpha);

bmw2008 = (manulist == 1);
volks2008 = (manulist == 21);
nbmw2008 = sum(bmw2008);
nvolks2008 = sum(volks2008);

a_pre_v_price = prices2008'*volks2008/nvolks2008;
a_pre_b_price = prices2008'*bmw2008/nbmw2008;
a_pre_v_mp = (pcd2008./C2008)'*volks2008/nvolks2008;
a_pre_b_mp = (pcd2008./C2008)'*bmw2008/nbmw2008;
t_pre_v_profit = profits2008'*volks2008; 
t_pre_b_profits = profits2008'*bmw2008;

premerger = [a_pre_price a_pre_cost a_pre_mp t_pre_ind a_pre_cs ...
             a_pre_v_price a_pre_b_price a_pre_v_mp a_pre_b_mp  ...
             t_pre_v_profit t_pre_b_profits]';


%% Diagnostics
% hist(C2008,100)
% hist(markup2008,100)
% 
%%
%% 1.1) Merger and Economies of Scale

ownershipmatrix_merger1 = zeros(length(market2008));


for i = 1:length(market2008)
    if ismember(i,find(manulist == 1)) || ismember(i,find(manulist == 21))
        ownershipmatrix_merger1(:,i) = (manulist == manulist(1)) + (manulist == manulist(964));
    else
        ownershipmatrix_merger1(:,i) = (manulist == manulist(i));
    end    
end

opt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 40000,'MaxIter',5000);

price_init = prices2008;
p_merger = fsolve(@(price)sys_foc(price, delta2008_np, alpha, C2008, ownershipmatrix_merger1, nu, theta2), price_init, opt);

pcd_merger = p_merger-C2008;
[~,shares_merger] = omegafinder(delta2008_np, theta2, p_merger, nu, alpha);
profits_merger = pcd_merger.*shares_merger*num_hh2008;

a_post_price = mean(p_merger);
a_post_cost = mean(C2008);
a_post_mp = mean(pcd_merger./(C2008));
t_post_ind = sum(profits_merger,1);
a_post_cs = consumer(delta2008_np,theta2,p_merger,nu,alpha);

a_post_v_price = p_merger'*volks2008/nvolks2008;
a_post_b_price = p_merger'*bmw2008/nbmw2008;
a_post_v_mp = (pcd_merger./C2008)'*volks2008/nvolks2008;
a_post_b_mp = (pcd_merger./C2008)'*bmw2008/nbmw2008;
t_post_v_profit = profits_merger'*volks2008; 
t_post_b_profits = profits_merger'*bmw2008;

postmerger = [a_post_price a_post_cost a_post_mp t_post_ind a_post_cs ...
             a_post_v_price a_post_b_price a_post_v_mp a_post_b_mp  ...
             t_post_v_profit t_post_b_profits]';

pre_post = [premerger, postmerger];

%% 1.2) Efficiency Gains



efficiency_gains = 0.01:0.001:0.02;
gridsize = size(efficiency_gains,2);
cost_reduction = (efficiency_gains'*(C2008.*merging_parties1)')';

cost_grid = repmat(C2008,1,gridsize)-cost_reduction;
consumer_surplus_grid = zeros(gridsize,1);

for costindex = 1:gridsize   
    p_merger = fsolve(@(price)sys_foc(price, delta2008_np, alpha, cost_grid(:,costindex), ownershipmatrix_merger1, nu, theta2), price_init, opt);
    consumer_surplus_grid(costindex) = consumer(delta2008_np, theta2, p_merger, nu, alpha);
end

[min_cs, loc_min] = min(consumer_surplus_grid);
eff_gains_min = efficiency_gains(loc_min);

%% 2.1) One Way Cross Participation

ownershipmatrix_RtoP = zeros(length(market2008));

for i = 1:length(market2008)
    if ismember(i,find(manulist == 15)) 
        ownershipmatrix_RtoP(i,:) = (manulist == manulist(i)) + 0.3*(manulist == 14);
    else
        ownershipmatrix_RtoP(i,:) = (manulist == manulist(i));
    end    
end
psa2008 = (manulist == 14);
renault2008 = (manulist == 15);

cr_base_aprice = mean(prices2008);
cr_base_amp = mean((prices2008-C2008)./C2008);
cr_base_profits = (prices2008-C2008).*shares2008*num_hh2008;
cr_base_r_profits = cr_base_profits'*renault2008;
cr_base_p_profits = cr_base_profits'*psa2008;
cr_base_r_profits_adj = cr_base_r_profits;
cr_base_p_profits_adj = cr_base_p_profits;
cr_base_ind_profits = sum(cr_base_profits,1);
cr_base_cs = consumer(delta2008_np, theta2, prices2008, nu, alpha);

baseline = [cr_base_aprice cr_base_amp cr_base_r_profits ...
            cr_base_p_profits cr_base_r_profits_adj cr_base_p_profits_adj ...
            cr_base_ind_profits cr_base_cs];

price_init = prices2008;
p_RtoP = fsolve(@(price)sys_foc(price, delta2008_np, alpha, C2008, ownershipmatrix_RtoP, nu, theta2), price_init, opt);
[~,shares_RtoP] = omegafinder(delta2008_np, theta2, p_RtoP, nu, alpha);

cr_RtoP_aprice = mean(p_RtoP);
cr_RtoP_amp = mean((p_RtoP-C2008)./C2008);
cr_RtoP_profits = (p_RtoP-C2008).*shares_RtoP*num_hh2008;
cr_RtoP_r_profits = cr_RtoP_profits'*renault2008;
cr_RtoP_p_profits = cr_RtoP_profits'*psa2008;
cr_RtoP_r_profits_adj = cr_RtoP_r_profits+0.3*cr_RtoP_p_profits;
cr_RtoP_p_profits_adj = 0.7*cr_RtoP_p_profits;
cr_RtoP_ind_profits = sum(cr_RtoP_profits,1);
cr_RtoP_cs = consumer(delta2008_np, theta2, p_RtoP, nu, alpha);

RtoP = [cr_RtoP_aprice cr_RtoP_amp cr_RtoP_r_profits ...
        cr_RtoP_p_profits cr_RtoP_r_profits_adj cr_RtoP_p_profits_adj ...
        cr_RtoP_ind_profits cr_RtoP_cs];


%% 2.2) Two way Cross Participation

%% PSA buys Renault

ownershipmatrix_PtoR = zeros(length(market2008));

for i = 1:length(market2008)
    if ismember(i,find(manulist == 15)) 
        ownershipmatrix_PtoR(i,:) = (manulist == manulist(i)) + 0.3*(manulist == 14);
    elseif ismember(i,find(manulist == 14))
        ownershipmatrix_PtoR(i,:) = (manulist == manulist(i)) + 0.3*(manulist == 15);
    else
        ownershipmatrix_PtoR(i,:) = (manulist == manulist(i));
    end    
end

price_init = prices2008;
p_PtoR = fsolve(@(price)sys_foc(price, delta2008_np, alpha, C2008, ownershipmatrix_PtoR, nu, theta2), price_init, opt);
[~,shares_PtoR] = omegafinder(delta2008_np, theta2, p_PtoR, nu, alpha);

cr_PtoR_aprice = mean(p_PtoR);
cr_PtoR_amp = mean((p_PtoR-C2008)./C2008);
cr_PtoR_profits = (p_PtoR-C2008).*shares_PtoR*num_hh2008;
cr_PtoR_r_profits = cr_PtoR_profits'*renault2008;
cr_PtoR_p_profits = cr_PtoR_profits'*psa2008;
cr_PtoR_r_profits_adj = 0.7*cr_PtoR_r_profits+0.3*cr_PtoR_p_profits;
cr_PtoR_p_profits_adj = 0.7*cr_PtoR_p_profits+0.3*cr_PtoR_r_profits;
cr_PtoR_ind_profits = sum(cr_PtoR_profits,1);
cr_PtoR_cs = consumer(delta2008_np, theta2, p_PtoR, nu, alpha);

PtoR = [cr_PtoR_aprice cr_PtoR_amp cr_PtoR_r_profits ...
        cr_PtoR_p_profits cr_PtoR_r_profits_adj cr_PtoR_p_profits_adj ...
        cr_PtoR_ind_profits cr_PtoR_cs];


%% PSA buys Renault before Renault buys PSA

ownershipmatrix_PtoRnoR = zeros(length(market2008));

for i = 1:length(market2008)
    if ismember(i,find(manulist == 14)) 
        ownershipmatrix_PtoRnoR(i,:) = (manulist == manulist(i)) + 0.3*(manulist == 15);
    else
        ownershipmatrix_PtoRnoR(i,:) = (manulist == manulist(i));
    end    
end

price_init = prices2008;
p_PtoRnoR = fsolve(@(price)sys_foc(price, delta2008_np, alpha, C2008, ownershipmatrix_PtoRnoR, nu, theta2), price_init, opt);
[~,shares_PtoRnoR] = omegafinder(delta2008_np, theta2, p_PtoRnoR, nu, alpha);

cr_PtoRnoR_aprice = mean(p_PtoRnoR);
cr_PtoRnoR_amp = mean((p_PtoRnoR-C2008)./C2008);
cr_PtoRnoR_profits = (p_PtoRnoR-C2008).*shares_PtoRnoR*num_hh2008;
cr_PtoRnoR_r_profits = cr_PtoRnoR_profits'*renault2008;
cr_PtoRnoR_p_profits = cr_PtoRnoR_profits'*psa2008;
cr_PtoRnoR_r_profits_adj = 0.7*cr_PtoRnoR_r_profits;
cr_PtoRnoR_p_profits_adj = cr_PtoRnoR_p_profits+0.3*cr_PtoRnoR_r_profits;
cr_PtoRnoR_ind_profits = sum(cr_PtoRnoR_profits,1);
cr_PtoRnoR_cs = consumer(delta2008_np, theta2, p_PtoRnoR, nu, alpha);

PtoRnoR = [cr_PtoRnoR_aprice cr_PtoRnoR_amp cr_PtoRnoR_r_profits ...
        cr_PtoRnoR_p_profits cr_PtoRnoR_r_profits_adj cr_PtoRnoR_p_profits_adj ...
        cr_PtoRnoR_ind_profits cr_PtoRnoR_cs];

cr = [baseline; RtoP; PtoR; PtoRnoR ]';


























