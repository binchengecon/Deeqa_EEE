rng(2021);


T = 1200;

% Set parameters
%%demand parameters
a0 = 20;
a1 = 4;
a2 = 1;
b0 = 1;
%b1 = b2 =0

%%marginal cost parameters
m0 = 2;
m1 = 2;
m2 = 1;
m3 = 1;

%%stack true parameters in a column vector
beta_d_true = [a0;a1;a2;b0];

beta_s_true = [m0; m1; m2; m3];


% Draw the exogenous variables

pop = rand(T,1);%Uniform(0,1)

inc =  4*rand(T,1);%Uniform(0,4)

cap =  10+10*rand(T,1);%Uniform(10,20)

wage = 2*rand(T,1);%Uniform(0,2)

material = 3*rand(T,1);%U(0,3)

% Draw the shocks from N(0,1)
eps_d = randn(T,1);%demand shock

eps_s = randn(T,1);%supply shock

% Fixed costs, represents the sum of firms' fixed costs
ft = 2+ 10*rand(T,1); %U(2,12) implies f0 = 7

f0 = 7; %parameter of fixed cost

%%Combine exogenous variables and parameters

%exogenous part of the marginal cost: mc = mc_tilde+m3*Qt/Kt
mc_tilde = m0 + m1*wage + m2*material + eps_s;

%%demand
at = a0 + a1*pop + a2*inc + eps_d;

bt = b0*ones(T,1);

%%Solve for model equilibrium: Nt, Pt, Qt

%initialize arrays
data = zeros(T,2); %to store prices (1st column) and quantities (2nd column)

n_firm = zeros(T,1); %to store the Nt the number of firms

profit = zeros(T,1); %store profits

mc = zeros(T,1); %store full marginal cost:  mc = mc_tilde+m3*Qt/Kt

%Loop over markets, for each market we need to solve for {Nt, Pt, Qt}
for i = 1:T
    
    profit_temp = 10;%initialize the loop
    
    n_temp = 0;

    %while loop to check if profits is positive for N=1,... firms. 
    %"while loop" stops whenever the condition "profit_temp >= 0" is not satisfied
    %need to intialize the value profit_temp, to a positive value to start
    %the loop
    %need to intialize the number of firms "n_temp" to 0, at each iteration
    %of the while loop, n_temp is incremented by 1
    
    while profit_temp >= 0 %condition of the while loop
        n_temp = n_temp+1; %increment number of firms
        
        %write the matrix with the parameters of the system of equation
        %defining the optimal prices and quantities
        %The system is A*(P;Q) = B
        A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/n_temp];
        B = [at(i); mc_tilde(i)];
        
        endog = A \ B; %solve the system of FOC in (P,Q) using a matrix 
        %inversion: A\B = A^-1*B, "\" is the "leftdivide" matrix function,
        %more efficient than the function "inv()"
        
        mc_temp = mc_tilde(i) + m3*endog(2)/cap(i); %calculate the full 
        %marginal cost given the number of firm
        
        %compute the total profit given the prices and quantities
        profit_temp = (endog(1) - mc_temp).*endog(2) - ft(i); %(Pt-mc)*Qt-ft
    end
    %since the while loop stops when we have negative profits, the
    %equilibrium number of firms is the previous number minus 1
    n_firm(i) = n_temp-1;
    
    %solve (again!) for equilibrium prices
    A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/(n_firm(i))];
    B = [at(i); mc_tilde(i)];
    endog = A \ B;
    %store equilibrium price and quantity
    data(i,:) = endog';
    %store full marginal cost
    mc(i) = mc_tilde(i) + m3*endog(2)/cap(i);
    %store profit
    profit_temp = (endog(1) - mc(i)).*endog(2) - ft(i);
    profit(i) = profit_temp;
    
end

figure
histogram(n_firm)
title('Distribution of number of firms across markets')
xlabel('Number of firms')


if min(n_firm) == 0
    disp('Number of firms = 0')
    pause
end


%%end of generating data

%% Estimation part: assume we do not know the parameters and want to estimate them

price = data(:,1); %Pt: the first column of data
qty = data(:,2); %Qt: the second column of data

p_check = at-bt.*qty; %check if price consistent with demand condition

if max(abs(p_check-price)) > 1e-4
    disp('Check prices')
    pause
end
    
%generate the matrix of exogenous regressors in the demand equation
xo_d = [ones(T,1), pop, inc];

%generate the matrix of endogenous regressors in the demand equation
xo_d_endo = [qty];

%generate the matrix of exogenous regressors in the supply equation
xo_s = [ones(T,1), wage, material];

%generate the matrix of endogenous regressors in the supply equation
xo_s_endo = qty./cap;


%%OLS estimation
beta_d_ols = regress(price, [xo_d, -xo_d_endo]);

%%IV estimation
%generate matrix of instruments for demand
z_d = [wage, cap, material];

%generate matrix of instruments for supply
z_s = [pop, inc];

[beta_d_iv, eps_hat_d] = ivregress(price, -xo_d_endo, xo_d, z_d);

bt_hat = beta_d_iv(4); %store bt_hat

%Estimate supply

%first calculate the dependent variable mc from estimates of demand
%parameter bt_hat

mc_hat = price - bt_hat*qty./n_firm ;

%OLS estimation
beta_s_ols = regress(mc_hat, [xo_s, xo_s_endo]);

%IV estimation
[beta_s_iv, eps_hat_s] = ivregress(mc_hat, xo_s_endo, xo_s, z_s);


%Store all estimates and compare with true parameters
%in order: true, OLS, IV
disp('compare demand parameters: True, OLS, IV')
[beta_d_true, beta_d_ols, beta_d_iv]

disp('compare supply parameters: True, OLS, IV')
[beta_s_true, beta_s_ols, beta_s_iv]


%estimate of f0: first, compute profits

%compute estimates of At, don't forget the residual which is part of the
%demand function (= demand shock) and Bt
at_hat = xo_d*beta_d_iv(1:3) + eps_hat_d;

bt_hat = beta_d_iv(4)*ones(T,1);

%get calcualte the marginal costs given the estimated parameters, including
%the cost shocks
m3_hat = beta_s_iv(4);

mc_tilde_hat = xo_s * beta_s_iv(1:3) + eps_hat_s;

%compute profits given marginal cost
profits_hat = (price - mc_hat).*qty;

%profit if there is an additional firm
%initialize vector of counterfactual profit
profits_n_plus_one = zeros(T,1);

%loop over market
for i = 1:T
    %solve for counterfactual price and quantity (if there is N+1 firms on
    %the market)
    A = [1, bt_hat(i) ; 1, -m3_hat/cap(i)- bt_hat(i)/(n_firm(i)+1)];
    B = [at_hat(i); mc_tilde_hat(i)];
    endog = A \ B;
    
    %calculate marginal cost and profit
    mc_temp = mc_tilde_hat(i) + m3_hat*endog(2)./cap(i);
    profits_n_plus_one(i) = (endog(1) - mc_temp).*endog(2);
end
    
%check that profit(N) > profit(N+1)
if min(profits_hat-profits_n_plus_one)<0
    disp('problem bounds')
    pause
end


%Estimate of fixed cost as the middle of the interval
f0_hat = mean((profits_hat + profits_n_plus_one)/2);

%compare estimate with true:
disp('compare true f0 with estimate')
[f0, f0_hat]

%% counterfactual simulations
%Predict the effect of a green tax of sales: price_paid = price_firm*(1.2)
%initialize arrays
data_counter = zeros(T,2);
n_firm_counter = zeros(T,1);
profit_counter = zeros(T,1);

%loop over markets
for i = 1:T
    %initialize profits and number of firms
    profit_temp = 10;
    n_temp = 0;
    %while profit is positive increase the number of firms
    while profit_temp >= 0
        n_temp = n_temp+1;
        %here coefficient of pt is multiplied by 1.2
        A = [1.2, bt_hat(i) ; 1, -m3_hat/cap(i)- bt_hat(i)/n_temp];
        B = [at_hat(i); mc_tilde_hat(i)];
        endog = A \ B;
        %calculate marginal costs and industry profit given the estimated
        %fixed costs
        mc_temp = mc_tilde_hat(i) + m3_hat*endog(2)/cap(i);
        profit_temp = (endog(1) - mc_temp).*endog(2) - f0_hat;
    end
    %once profit becomes negative take the previous number of firms and
    %recalculate prices and quantities
    n_firm_counter(i) = n_temp-1;
    A = [1, bt_hat(i) ; 1, -m3_hat/cap(i)- bt_hat(i)/(n_firm_counter(i))];
    B = [at_hat(i); mc_tilde_hat(i)];
    endog = A \ B;
    data_counter(i,:) = endog';
    mc_temp = mc_tilde_hat(i) + m3_hat*endog(2)/cap(i);
    profit_temp = (endog(1) - mc_temp).*endog(2) - f0_hat;
    profit_counter(i) = profit_temp;
    
end

figure
histogram(n_firm_counter)
title('Counterfactual distribution of number of firms across markets')
xlabel('Number of firms')
axis([0 100 0 125])

if min(n_firm_counter) == 0
    disp('check number of firms in counterfactual')
    pause
end

%% To make a fair comparison, need to compare the counterfactual number of firms with the 
%predicted number of firms using our model

%initialize arrays
data_counter_initial = zeros(T,2);
n_firm_counter_initial = zeros(T,1);
profit_counter_initial = zeros(T,1);

%loop over markets
for i = 1:T
    %initialize profits and number of firms
    profit_temp = 10;
    n_temp = 0;
    %while profit is positive increase the number of firms
    while profit_temp >= 0
        n_temp = n_temp+1;
        %write system that defines equilibrium price and quantity
        A = [1, bt_hat(i) ; 1, -m3_hat/cap(i)- bt_hat(i)/n_temp];
        B = [at_hat(i); mc_tilde_hat(i)];
        endog = A \ B;
        %calculate marginal costs and industry profit given the estimated
        %fixed costs
        mc_temp = mc_tilde_hat(i) + m3_hat*endog(2)/cap(i);
        profit_temp = (endog(1) - mc_temp).*endog(2) - f0_hat;
    end
    %once profit becomes negative take the previous number of firms and
    %recalculate prices and quantities
    n_firm_counter_initial(i) = n_temp-1;
    A = [1, bt_hat(i) ; 1, -m3_hat/cap(i)- bt_hat(i)/(n_firm_counter_initial(i))];
    B = [at_hat(i); mc_tilde_hat(i)];
    endog = A \ B;
    data_counter_initial(i,:) = endog';
    mc_temp = mc_tilde_hat(i) + m3_hat*endog(2)/cap(i);
    profit_temp = (endog(1) - mc_temp).*endog(2) - f0_hat;
    profit_counter_initial(i) = profit_temp;
    
end

figure
histogram(n_firm_counter_initial)
title('Counterfactual INITIAL distribution of number of firms across markets')
xlabel('Number of firms')
axis([0 100 0 125])

if min(n_firm_counter_initial) == 0
    disp('check number of firms in counterfactual initial')
    pause
end



%average number of firms before and after the tax
disp('average number of firms observed and predicted before and after the tax')
[mean(n_firm), mean(n_firm_counter_initial), mean(n_firm_counter)]


mean(n_firm_counter_initial < n_firm_counter)
%n_firm counter should be compared with predicted number of firms from the
%model. No market has a larger number of firms

%plot the two distributions of number of firms on the same graph: use
%kernel density estimation to smooth the histograms

[density_firms_init, x_axis_init] = ksdensity(n_firm_counter_initial);

[density_firms_tax, x_axis_counter] = ksdensity(n_firm_counter);

figure
plot(x_axis_init, density_firms_init, 'linewidth', 2)
hold on
plot(x_axis_counter,density_firms_tax, 'linewidth', 2)
hold off
legend('Predicted initial number of firms', 'Predicted with the tax')
xlabel('Number of firms')
ylabel('Density')

%% This is the same data generation but xi_d and xi_s have lower variances

% New shocks, draw them with lower variance
eps_d_lower = 0.5*randn(T,1);
eps_s_lower = 0.5*randn(T,1);

%%generate variable cost

mc_tilde_lower = m0 + m1*wage + m2*material + eps_s_lower;

%%demand
at_lower = a0 + a1*pop + a2*inc + eps_d_lower;

%solve for Qt, Pt, Nt
data_lower = zeros(T,2);
n_firm_lower = zeros(T,1);
profit_lower = zeros(T,1);
for i = 1:T
    profit_temp = 10;
    n_temp = 0;
    while profit_temp >= 0
        n_temp = n_temp+1;
        A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/n_temp];
        B = [at_lower(i); mc_tilde_lower(i)];
        endog = A \ B;
        mc_temp = mc_tilde_lower(i) + m3*endog(2)/cap(i);
        profit_temp = (endog(1) - mc_temp).*endog(2) - ft(i);
    end
    n_firm_lower(i) = n_temp - 1;
    A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/n_firm_lower(i)];
    B = [at_lower(i); mc_tilde_lower(i)];
    endog = A \ B;
    mc_temp = mc_tilde_lower(i) + m3*endog(2)/cap(i);
    profit_temp = (endog(1) - mc_temp).*endog(2) - ft(i);
    profit_lower(i) = profit_temp;
    data_lower(i,:) = endog';   
    
end


if min(n_firm_lower) == 0
    disp('number of firms equal 0 with lower variance')
    pause
end

price_lower = data_lower(:,1);
qty_lower = data_lower(:,2);

xo_d_endo_lower = [qty_lower];

xo_s_endo_lower = qty_lower./cap;


%%OLS
beta_d_ols_lower = regress(price_lower, [xo_d, -xo_d_endo]);

%%IV

[beta_d_iv_lower, eps_hat_d_lower] = ivregress(price_lower, -xo_d_endo_lower, xo_d, z_d);

%true, OLS, IV

bt_hat_lower = beta_d_iv_lower(4);
%Estimate supply
mc_hat_lower = price_lower - bt_hat_lower*qty_lower./n_firm_lower ;

beta_s_ols_lower = regress(mc_hat_lower, [xo_s, xo_s_endo_lower]);


[beta_s_iv_lower, eps_hat_s_lower] = ivregress(mc_hat_lower, xo_s_endo_lower, xo_s, z_s);

%comparison of estimates: true, OLS, IV
disp('compare demand estimates: true, OLS, IV')
beta_d_lowvar = [beta_d_true, beta_d_ols_lower, beta_d_iv_lower]

disp('compare supply estimates: true, OLS, IV')
beta_s_lowvar = [beta_s_true, beta_s_ols_lower, beta_s_iv_lower]

%estimate of f0
at_hat_lower = xo_d*beta_d_iv_lower(1:3) + eps_hat_d_lower;

bt_hat_lower = beta_d_iv_lower(4);

m3_hat_lower = beta_s_iv_lower(4);

mc_tilde_hat_lower = xo_s * beta_s_iv_lower(1:3) + eps_hat_s_lower;

profits_hat_lower = (price_lower - mc_hat_lower).*qty_lower;

%calculate variable part of the profit(N+1)
profits_n_plus_one_lower = zeros(T,1);
for i = 1:T
    %solve for Pt, Qt when we have N+1 firms
    A = [1, bt_hat_lower ; 1, -m3_hat_lower/cap(i) - bt_hat_lower/(n_firm_lower(i)+1)];
    B = [at_hat_lower(i); mc_tilde_hat_lower(i)];
    endog = A \ B;
    %compute full marginal cost
    mc_temp = mc_tilde_hat_lower(i) + m3_hat_lower*endog(2)./cap(i);
    profits_n_plus_one_lower(i) = (endog(1) - mc_temp).*endog(2);
end
    
%check that profit(N) > profit(N+1)
if min(profits_hat-profits_n_plus_one) < 0
    disp('problem bounds lower variance')
    pause
end


f0_hat_lowvar = mean((profits_hat_lower + profits_n_plus_one_lower)/2);


disp('compare demand estimates with lower variance')
%demand
[beta_d_true, beta_d_iv, beta_d_iv_lower]

disp('compare supply estimates with lower variance')
%supply
[beta_s_true, beta_s_iv, beta_s_iv_lower]

disp('compare fixed cost estimates with lower variance')
%fixed cost
[f0, f0_hat, f0_hat_lowvar]