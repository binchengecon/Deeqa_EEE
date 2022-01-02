clear all
cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Homework\Correction'
data = csvread('data2020.csv');

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

num_brand = data(:,3);

price = data(:,4);

xo = data(:,5:8);
%in order: cylinder, weight, HP, fuel cost

nb_x = size(xo,2) + 1;

no_hh = data(:,9);

sj = data(:,10);

iv = data(:,11:22);

s0_temp = 1-accumarray(num_market, sj);

s0 = s0_temp(num_market);

ye_fe = zeros(kt, T-1);
for tt = 2:T
    ye_fe(:,tt-1) = (num_market==tt);
end

brd_fe = zeros(kt, max(num_brand)-1);
for ii = 2: max(num_brand)
    brd_fe(:,ii-1) = (num_brand==ii);
end

y = log(sj./s0);


opt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');


rng(2020);
ns = 100;
nu = randn(1, ns);

%% Q1

xo_nofe = [ones(kt,1), xo];

z_nofe = [xo_nofe, iv];

n_z = size(z_nofe,2);

invZ_nofe = (z_nofe'*z_nofe)\eye(n_z);


[sigma_hat_nofe] = fminunc(@(theta)gmm_random(theta, xo_nofe, price, price, nu,...
    sj, z_nofe, invZ_nofe, y, marketStarts, marketEnds), -1.5, opt);

mu_est_nofe = sigma_hat_nofe*price*nu;

delta_est_nofe = zeros(kt,1);
for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    delta_est_nofe(Jt) = inv_sh_by_market(sj(Jt), y(Jt), mu_est_nofe(Jt,:));
end


[beta_hat_nofe] = ivregress(delta_est_nofe, xo_nofe, price, z_nofe, invZ_nofe);

alpha_hat_nofe = beta_hat_nofe(end);

parms_nofe = [alpha_hat_nofe;sigma_hat_nofe;beta_hat_nofe(1:nb_x)];


alpha_i_nofe = alpha_hat_nofe + nu*sigma_hat_nofe;

%plot the distribution of the price coefficient
[f_alpha_nofe, x_alpha_nofe]= ksdensity(alpha_i_nofe);
figure
plot(x_alpha_nofe, f_alpha_nofe)


%% Q2
xo_big = [ones(kt,1), xo, ye_fe, brd_fe];

z = [xo_big, iv];

n_z = size(z,2);

invZ = (z'*z)\eye(n_z);


grid = [-4:0.5:4];
nb_grid = size(grid,2);
fval_grid = zeros(nb_grid,1);
for ii = 1:nb_grid
    fval_grid(ii) = gmm_random(grid(ii), xo_big, price, price, nu,...
        sj, z, invZ, y, marketStarts, marketEnds);
    
end

%check the shape of the objective function
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

parms_fe = [alpha_hat;sigma_hat;beta_hat(1:nb_x)];
delta_hp = delta_est - alpha_hat*price;

alpha_i = alpha_hat + nu*sigma_hat;

%plot the distribution of the price coefficient
[f_alpha, x_alpha]= ksdensity(alpha_i);
figure
plot(x_alpha, f_alpha)

%compare with and without FE
figure
plot(x_alpha, f_alpha)
hold on
plot(x_alpha_nofe, f_alpha_nofe)
hold off


%get marginal costs
cost = zeros(kt,1);

for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    omeg = omega_rc(delta_est(Jt), price(Jt), num_group(Jt), alpha_hat, sigma_hat, nu);
    cost(Jt) = price(Jt) + omeg\sj(Jt);
end

%check costs are positive
if min(cost)<0
    disp('some costs negative')
else
    disp('all costs positive')
end



%check that we get same price from marginal costs
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


%% Q3 diff IV
nb_xo = size(xo,2);

std_x = zeros(T,nb_xo);
for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    std_x(tt,:) = std(xo(Jt,:));
end
iv_diff1 = zeros(kt,nb_xo);
iv_diff2 = zeros(kt,nb_xo);

for xx = 1:nb_xo
    for tt = 1:T
        Jt = marketStarts(tt):marketEnds(tt);
        
        for ii = Jt
            iv_diff1(ii, xx) = sum(abs(xo(ii, xx) - xo(Jt, xx)).* (abs(xo(ii,xx) - xo(Jt, xx)) <=std_x(tt, xx))) ;
            iv_diff2(ii,xx) = sum((xo(ii, xx) - xo(Jt, xx)).^2);
        end
    end
end

p_hat = z*((z'*z)\(z'*price));

phat_diff1 = zeros(kt, 1);
phat_diff2 = zeros(kt, 1);
for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    for ii = Jt
        phat_diff1(ii) = sum(abs(p_hat(ii) - p_hat(Jt)).* (abs(p_hat(ii) - p_hat(Jt)) <=std(p_hat(Jt)))) ;
        phat_diff2(ii) = sum((p_hat(ii) - p_hat(Jt)).^2)/1e5;
    end
end

z_new1 = [z, iv_diff1, phat_diff1];

z_new2 = [z, iv_diff2, phat_diff2];


invZ_new1 = (z_new1'*z_new1)\eye(size(z_new1, 2), size(z_new1,2));



grid = [-4:0.5:4];
nb_grid = size(grid,2);
fval_grid = zeros(nb_grid,1);
for ii = 1:nb_grid
    fval_grid(ii) = gmm_random(grid(ii), xo_big, price, price, nu,...
        sj, z_new1, invZ_new1, y, marketStarts, marketEnds);
    
end

figure
plot(grid, fval_grid)


[sigma_hat_new1, fval_new1] = fminunc(@(theta)gmm_random(theta, xo_big, price, price, nu,...
    sj, z_new1, invZ_new1, y, marketStarts, marketEnds), 1.5, opt);

mu_est = sigma_hat_new1*price*nu;

delta_est_new1 = zeros(kt,1);
for tt = 1:T
    Jt = marketStarts(tt):marketEnds(tt);
    delta_est_new1(Jt) = inv_sh_by_market(sj(Jt), y(Jt), mu_est(Jt,:));
end

[beta_hat_new1] = ivregress(delta_est_new1, xo_big, price, z_new1, invZ_new1);

alpha_hat_new1 = beta_hat_new1(end);

parms_new1 = [alpha_hat_new1; sigma_hat_new1; beta_hat_new1(1:nb_x)];

alpha_i = alpha_hat_new1 + nu*sigma_hat_new1;

%plot the distribution of the price coefficient
[f_alpha_new1, x_alpha_new1]= ksdensity(alpha_i);
figure
plot(x_alpha_new1, f_alpha_new1)

%compare with and without FE
figure
plot(x_alpha, f_alpha, 'lineWidth', 1.7)
hold on
plot(x_alpha_nofe, f_alpha_nofe, 'lineWidth', 1.7)
plot(x_alpha_new1, f_alpha_new1, 'lineWidth', 1.7)
xlabel('\alpha')
ylabel('Density')
legend('No FE', 'FE', 'G & H IV')
axis tight
hold off

% saveas(gcf, '/Users/isisdurrmeyer/Dropbox/Teaching_DEEQA_2020/Projects/Homework_2020/distri_alpha.pdf')
%% Part 1: Question 3.5

% Collect data into tables for exportation into Latex


Table_Regression = zeros(7,2);
Table_CFResult   = zeros(6,2);
% No FE case
% beta_hat_Nfe: 1, cy,w,hp,    fuelcost, price

Table_Regression(1:5,1) = beta_hat_nofe(1:5);
Table_Regression(6,1)   = sigma_hat_nofe;
Table_Regression(7,1)   = beta_hat_nofe(end);

%FE case
% beta_hat_fe: 1, cy,w,hp, fuelcost,   brand(40 brands, 39 brand_fe) , price
%              1  2-4         5          6-44,                          45

Table_Regression(1:5,2) = beta_hat(1:5);
Table_Regression(6,2) = sigma_hat;
Table_Regression(7,2)   = beta_hat(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initial situation
Jt = marketStarts(T):marketEnds(T);

%lerner and mark up rate for 2008
av_price = mean(price(Jt))*10000;
av_mark1 = mean((price(Jt)-cost(Jt))./price(Jt))*100;
av_mark2 = (mean(price(Jt)./cost(Jt))-1)*100;

%profits
profit_j =  (price(Jt)-cost(Jt)).*sj(Jt).*no_hh(Jt)*10000/1e6; %in million euros

profit_group = accumarray(num_group(Jt), profit_j);
%REN (15), psa (14), bmw (1), VW (21)
profit_ren = profit_group(15);
profit_psa = profit_group(14);
profit_bmw = profit_group(1);
profit_vw = profit_group(21);

profit_ren_psa = profit_ren + profit_psa;
profit_vw_bmw = profit_bmw + profit_vw;

tot_profit = sum(profit_group);

av_cs = surplus_by_market(delta_est(Jt), price(Jt), sigma_hat, nu, alpha_hat);

res_init = [av_price; av_mark1; av_mark2; profit_vw_bmw; profit_ren_psa; tot_profit; av_cs];

%% Q1: merger between VW (21) & BMW (1)

num_group_bis = num_group(Jt);

num_group_bis(num_group(Jt) == 21) = 1;

sum(num_group_bis==21)

price_new1 = fsolve(@(price_temp) foc_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    num_group_bis, alpha_hat, sigma_hat, nu), price(Jt), opt);

delta_new1 = delta_hp(Jt) + alpha_hat*price_new1;

mu_new1 = sigma_hat*price_new1*nu;

sj_new1 = agg_sh_by_market(delta_new1, mu_new1);

av_price_new1 = mean(price_new1)*10000;

%post-merger Lerner index and mark up rate
av_mark1_new1 = mean((price_new1-cost(Jt))./price_new1)*100;
av_mark2_new1 = (mean(price_new1./cost(Jt))-1)*100;

profit_j_new1 =  (price_new1-cost(Jt)).*sj_new1.*no_hh(Jt)*10000/1e6;

profit_group_new1 = accumarray(num_group_bis, profit_j_new1);
%REN (15), psa (14), bmw-VW (1)

profit_vw_bmw_new1 = profit_group_new1(1);

tot_profit_new1 = sum(profit_group_new1);

av_cs_new1 = surplus_by_market(delta_new1, price_new1, sigma_hat, nu, alpha_hat);


res_new1 = [av_price_new1; av_mark1_new1; av_mark2_new1; profit_vw_bmw_new1; 0 ; tot_profit_new1; av_cs_new1];


%% Q2: min efficiency to have DCS = 0

%grid from 1% to 5% cost efficiency
grid = [0.01:0.005:0.05];
nb_grid = size(grid, 2);
csurplus = zeros(nb_grid,1);

for ii = 1:nb_grid
    
    cost_new = cost(Jt);
    
    cost_new(num_group_bis == 1) = (1-grid(ii)) * cost_new(num_group_bis == 1);
    
    price_new11 = fsolve(@(price_temp) foc_rc(price_temp, cost_new, delta_hp(Jt), ...
        num_group_bis, alpha_hat, sigma_hat, nu), price(Jt), opt);
    
    delta_new11 = delta_hp(Jt)+alpha_hat*price_new11;
    
    csurplus(ii) = surplus_by_market(delta_new11, price_new11, sigma_hat, nu, alpha_hat);
    
end

min_grid_refine = grid(max(find(csurplus<=av_cs)));

max_grid_refine = grid(min(find(csurplus>=av_cs)));

grid1 = grid*100;

figure
plot([min(grid1); max(grid1)], [av_cs; av_cs], 'linewidth', 1.7)
hold on
plot(grid1, csurplus, 'linewidth', 1.7)
xlabel('Efficiency gains (in %)')
ylabel('Consumer surplus')
axis tight
legend('Initial CS', 'location', 'best')
hold off

saveas(gcf, '/Users/isisdurrmeyer/Dropbox/Teaching_DEEQA_2020/Projects/Homework_2020/eff_gains.pdf')


%refining the grid
grid = [min_grid_refine:0.001:max_grid_refine];
nb_grid = size(grid, 2);
csurplus = zeros(nb_grid,1);

for ii = 1:nb_grid
    
    cost_new = cost(Jt);
    
    cost_new(num_group_bis == 1) = (1-grid(ii)) * cost_new(num_group_bis == 1);
    
    price_new11 = fsolve(@(price_temp) foc_rc(price_temp, cost_new, delta_hp(Jt), ...
        num_group_bis, alpha_hat, sigma_hat, nu), price(Jt), opt);
    
    delta_new11 = delta_hp(Jt)+alpha_hat*price_new11;
    
    csurplus(ii) = surplus_by_market(delta_new11, price_new11, sigma_hat, nu, alpha_hat);
    
end

min_grid_refine2 = grid(max(find(csurplus<=av_cs)));

max_grid_refine2 = grid(min(find(csurplus>=av_cs)));

%efficiency between min & max, take average

eff_gains = (min_grid_refine2 + max_grid_refine2)/2





%% Q2 effect of cross participations 
%Q2-1) REN (15) buys 30% PSA (14)
%i.e. REN internalizes 30% of PSA's profits
kkt = size(Jt,2);
kappa_mat = zeros(kkt,kkt);

kappa_mat(num_group(Jt)==15, num_group(Jt)==14) = 0.3;

mat_own = repmat(num_group(Jt),1,kkt)==repmat(num_group(Jt)',kkt,1);

owner = mat_own + kappa_mat;

price_new4 = fsolve(@(price_temp) foc_kappa_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    alpha_hat, sigma_hat, nu, owner), price(Jt), opt);


delta_new4 = delta_hp(Jt) + alpha_hat*price_new4;

mu_new4 = sigma_hat*price_new4*nu;

[sj_new4] = agg_sh_by_market(delta_new4, mu_new4);

%post-merger price, Lerner index and mark up rate
av_price_new4 = mean(price_new4)*10000;
av_mark1_new4 = mean((price_new4-cost(Jt))./price_new4)*100;
av_mark2_new4 = (mean(price_new4./cost(Jt))-1)*100;


%profits in 2008
profit_j_new4 =  (price_new4-cost(Jt)).*sj_new4.*no_hh(Jt)*10000/1e6;

profit_group = accumarray(num_group(Jt), profit_j_new4);

%REN (15), psa (14), bmw (1), VW (21)
profit_ren_new4 = profit_group(15);
profit_psa_new4 = profit_group(14);

%true profit, accounting for participation
profit_ren_alt_new4 = profit_group(15)+0.3*profit_group(14);
profit_psa_alt_new4 = 0.7* profit_group(14);

profit_ren_psa_new4 = profit_ren_new4 + profit_psa_new4;

tot_profit_new4 = sum(profit_group);

av_cs_new4 = surplus_by_market(delta_new4, price_new4, sigma_hat, nu, alpha_hat);

res_new4 = [av_price_new4; av_mark1_new4; av_mark2_new4; 0; profit_ren_psa_new4;...
    tot_profit_new4; av_cs_new4];



%% Q2-2) PSA buys 30% REN in exchange
kappa_mat(num_group(Jt)==14, num_group(Jt)==15) = 0.3/0.7;
kappa_mat(num_group(Jt)==15, num_group(Jt)==14) = 0.3/0.7;

owner22 = mat_own + kappa_mat;

price_new5 = fsolve(@(price_temp) foc_kappa_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    alpha_hat, sigma_hat, nu, owner22), price_new4, opt);


delta_new5 = delta_hp(Jt) + alpha_hat*price_new5;

mu_new5 = sigma_hat*price_new5*nu;

[sj_new5] = agg_sh_by_market(delta_new5, mu_new5);

%post-merger price, Lerner index and mark up rate
av_price_new5 = mean(price_new5)*10000;
av_mark1_new5 = mean((price_new5-cost(Jt))./price_new5)*100;
av_mark2_new5 = (mean(price_new5./cost(Jt))-1)*100;


%profits in 2008
profit_j_new5 =  (price_new5-cost(Jt)).*sj_new5.*no_hh(Jt)*10000/1e6;

profit_group = accumarray(num_group(Jt), profit_j_new5);

%REN (15), psa (14), bmw (1), VW (21)
profit_ren_new5 = profit_group(15);
profit_psa_new5 = profit_group(14);

%true profit, accounting for participation
profit_ren_alt_new5 = 0.7*profit_group(15)+0.3*profit_group(14);
profit_psa_alt_new5 = 0.7* profit_group(14)+0.3*profit_group(15);

profit_ren_psa_new5 = profit_ren_new5 + profit_psa_new5;

tot_profit_new5 = sum(profit_group);

av_cs_new5 = surplus_by_market(delta_new5, price_new5, sigma_hat, nu, alpha_hat);

res_new5 = [av_price_new5; av_mark1_new5; av_mark2_new5; 0; profit_ren_psa_new5;...
    tot_profit_new5; av_cs_new5];


%% Q2-2- PSA buys 30% REN , REN does not buy PSA
kappa_mat = zeros(kkt,kkt);
kappa_mat(num_group(Jt)==14, num_group(Jt)==15) = 0.3;

owner21 = mat_own + kappa_mat;

price_new5b = fsolve(@(price_temp) foc_kappa_rc(price_temp, cost(Jt), delta_hp(Jt), ...
    alpha_hat, sigma_hat, nu, owner21), price_new4, opt);

delta_new5b = delta_hp(Jt) + alpha_hat*price_new5b;

mu_new5b = sigma_hat*price_new5b*nu;

[sj_new5b] = agg_sh_by_market(delta_new5b, mu_new5b);

%post-merger price, Lerner index and mark up rate
av_price_new5b = mean(price_new5b)*10000;
av_mark1_new5b = mean((price_new5b-cost(Jt))./price_new5b)*100;
av_mark2_new5b = (mean(price_new5b./cost(Jt))-1)*100;


%profits in 2008
profit_j_new5b =  (price_new5b-cost(Jt)).*sj_new5b.*no_hh(Jt)*10000/1e6;

profit_group = accumarray(num_group(Jt), profit_j_new5b);

%REN (15), psa (14), bmw (1), VW (21)
profit_ren_new5b = profit_group(15);
profit_psa_new5b = profit_group(14);

%true profit, accounting for participation
profit_ren_alt_new5b = 0.7*profit_group(15);
profit_psa_alt_new5b = profit_group(14)+0.3*profit_group(15);

profit_ren_psa_new5b = profit_ren_new5b + profit_psa_new5b;

tot_profit_new5b = sum(profit_group);

av_cs_new5b = surplus_by_market(delta_new5b, price_new5b, sigma_hat, nu, alpha_hat);

res_new5b = [av_price_new5b; av_mark1_new5b; av_mark2_new5b; 0; profit_ren_psa_new5b;...
    tot_profit_new5b; av_cs_new5b];


%% compute WTP for 30% PSA
matrix_profit_psa = [[profit_psa, profit_psa_alt_new4];[profit_psa_alt_new5b, profit_psa_alt_new5]];
%in order:
%No one buys, REN buys
%   PSA buys, both buy
