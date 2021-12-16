clearvars
clc

rng(2021);

T = 1200;
index = (1:1200)';

% Parameters
%%%demand parameters
a0 = 20;
a1 = 4;
a2 = 1;
b0 = 1;

%%cmarginal cost parameters
m0 = 2;
m1 = 2;
m2 = 1;
m3 = 1;

beta_s_true = [m0; m1; m2; m3];
beta_d_true = [a0;a1;a2;b0];


%% evaluate the bias and precision through MC: 100 sythetic datasets
ns = 30;
beta_s_ols_simu = zeros(4, ns);
beta_s_iv_simu = zeros(4, ns);
beta_d_ols_simu = zeros(4, ns);
beta_d_iv_simu = zeros(4, ns);
f0_hat_simu = zeros(ns,1);

for jj = 1:ns
    jj
    wage = 2*rand(T,1);
    
    pop = rand(T,1);
    
    cap =  10+10*rand(T,1);
    
    inc =  4*rand(T,1);
    
    input = 3*rand(T,1);
    
    eps_d = randn(T,1);
    eps_s = randn(T,1);
    
    ft = 2+ 10*rand(T,1); %U(2,12)
    f0 = 7;
    %%generate  variable cost
    
    mc_tilde = m0 + m1*wage + m2*input + eps_s;
    
    %% demand
    at = a0 + a1*pop + a2*inc + eps_d;
    
    bt = b0*ones(T,1);
    
    
    data = zeros(T,2);
    n_firm = zeros(T,1);
    profit = zeros(T,1);
    for i = 1:T
        profit_temp = 10;
        n_temp = 0;
        while profit_temp >= 0
            n_temp = n_temp+1;
            A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/n_temp];
            B = [at(i); mc_tilde(i)];
            endog = A \ B;
            mc_temp = mc_tilde(i) + m3*endog(2)/cap(i);
            profit_temp = (endog(1) - mc_temp).*endog(2) - ft(i);
        end
        n_firm(i) = n_temp-1;
        A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/(n_firm(i))];
        B = [at(i); mc_tilde(i)];
        endog = A \ B;
        data(i,:) = endog';
        mc_temp = mc_tilde(i) + m3*endog(2)/cap(i);
        profit_temp = (endog(1) - mc_temp).*endog(2) - ft(i);
        profit(i) = profit_temp;
        
    end
    
    if min(n_firm) == 0
        disp('check number of firms')
        pause
    end
    
    price = data(:,1);
    qty = data(:,2);
    
    p_check = at-bt.*data(:,2);
    
    xo_d = [ones(T,1), pop, inc];
    xo_d_endo = [qty];
    
    xo_s = [ones(T,1), wage, input];
    
    xo_s_endo = qty./cap;
    
    
    %%OLS
    beta_d_ols_simu(:,jj) = regress(price, [xo_d, -xo_d_endo]);
    
    %%IV
    z_d = [wage, cap, input];
    z_s = [pop, inc];
    
    [beta_d_iv_simu(:,jj), eps_hat_d] = ivregress(price, -xo_d_endo, xo_d, z_d);
    
    %true, OLS, IV
    
    bt_hat = beta_d_iv_simu(4,jj);
    %Estimate supply
    mc_hat = price - bt_hat*qty./n_firm ;
    
    beta_s_ols_simu(:,jj) = regress(mc_hat, [xo_s, xo_s_endo]);
    
    
    [beta_s_iv_simu(:,jj), eps_hat_s] = ivregress(mc_hat, xo_s_endo, xo_s, z_s);
    
    %estimate of f0: first, compute profits

    at_hat = xo_d*beta_d_iv_simu(1:3,jj) + eps_hat_d;

    bt_hat = beta_d_iv_simu(4,jj)*ones(T,1);

    m3_hat = beta_s_iv_simu(4,jj);

    mc_tilde_hat = xo_s * beta_s_iv_simu(1:3,jj) + eps_hat_s;

    profits_hat = (price - mc_hat).*qty;


    profits_n_plus_one = zeros(T,1);
    for i = 1:T
    A = [1, bt_hat(i) ; 1, -m3_hat/cap(i)- bt_hat(i)/(n_firm(i)+1)];
    B = [at_hat(i); mc_tilde_hat(i)];
    endog = A \ B;
    mc_temp = mc_tilde_hat(i) + m3_hat*endog(2)./cap(i);
    profits_n_plus_one(i) = (endog(1) - mc_temp).*endog(2);
    end
    

    f0_hat_simu(jj) = mean((profits_hat + profits_n_plus_one)/2);


    
    
    
end

%mean estimates and standard deviation
res_demand_simu = [beta_d_true, mean(beta_d_iv_simu,2), std(beta_d_iv_simu,[],2)]

res_supply_simu = [beta_s_true, mean(beta_s_iv_simu,2), std(beta_s_iv_simu,[],2)]

res_fixed_cost = [f0, mean(f0_hat_simu), std(f0_hat_simu)]



