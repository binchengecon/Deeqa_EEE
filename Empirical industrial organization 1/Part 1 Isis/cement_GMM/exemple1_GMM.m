rng(2020);


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

beta_d_true = [a0;a1;a2;b0];

beta_s_true = [m0; m1; m2; m3];


% Exogenous variables
wage = 2*rand(T,1);%Uniform(0,2)

pop = rand(T,1);%Uniform(0,1)

cap =  10+10*rand(T,1);%Uniform(10,20)

inc =  4*rand(T,1);

input = 3*rand(T,1);

eps_d = randn(T,1);%N(0,1)
eps_s = randn(T,1);


ft = 2+ 10*rand(T,1); %U(2,12) implies f0=7
f0 = 7;

%%generate  variable cost

mc_tilde = m0 + m1*wage + m2*input + eps_s;

%%demand
at = a0 + a1*pop + a2*inc + eps_d;

bt = b0*ones(T,1);

data = zeros(T,2);
n_firm = zeros(T,1);
profit = zeros(T,1);
mc = zeros(T,1);
for i = 1:T
    profit_temp = 10;
    n_temp = 0;
    while profit_temp >= 0
        n_temp = n_temp+1;
        A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/n_temp];
        B = [at(i); mc_tilde(i)];
        endog = A \ B; %A^-1*B
        mc_temp = mc_tilde(i) + m3*endog(2)/cap(i);
        profit_temp = (endog(1) - mc_temp).*endog(2) - ft(i);
    end
    n_firm(i) = n_temp-1;
    A = [1, bt(i) ; 1, -m3/cap(i)- bt(i)/(n_firm(i))];
    B = [at(i); mc_tilde(i)];
    endog = A \ B;
    data(i,:) = endog';
    mc(i) = mc_tilde(i) + m3*endog(2)/cap(i);
    profit_temp = (endog(1) - mc(i)).*endog(2) - ft(i);
    profit(i) = profit_temp;
    
end

figure
hist(n_firm)

if min(n_firm) == 0
    disp('check no of firms =0')
    pause
end
%check



price = data(:,1); %take the first column of matrix data
qty = data(:,2);

p_check = at-bt.*qty;

mean(abs(p_check-price))

xo_d = [ones(T,1), pop, inc];
xo_d_endo = [qty];

xo_s = [ones(T,1), wage, input];

xo_s_endo = qty./cap;

%%IV
z_d = [wage, cap, input];
z_s = [pop, inc];

z = [xo_d, z_d];

[beta_d_iv, eps_hat_d] = ivregress2(price, -xo_d_endo, xo_d, z_d);

bt_hat = beta_d_iv(4);
%Estimate supply

mc_hat = price - bt_hat*qty./n_firm ;


[beta_s_iv, eps_hat_s] = ivregress2(mc_hat, xo_s_endo, xo_s, z_s);

%% Estimation with GMM

x_d = [ones(T,1), pop, inc, - qty];

x_s = [ones(T,1), wage, input, qty./cap];

w_tilde = (z'*z)\eye(size(z,2));

w = blkdiag(w_tilde,w_tilde);

theta_true = [beta_d_true;beta_s_true];

gmm_true = gmm_obj2(theta_true, price, qty, x_d, x_s, z, w, n_firm);


theta_init = theta_true+rand(size(theta_true,1),1);

opt = optimoptions('fminunc', 'Display','iter', 'tolX', 1e-6, 'tolFun', 1e-6, 'OptimalityTolerance', 1e-8);

theta_gmm = fminunc(@(theta)gmm_obj2(theta, price, qty, x_d, x_s, z, w, n_firm), theta_init, opt);

%finale estimates, in order: true, IV, GMM

beta_d_gmm = theta_gmm(1:size(x_d,2));

beta_s_gmm = theta_gmm(size(x_d,2)+1:end);

beta_d_comp = [beta_d_true, beta_d_iv, beta_d_gmm]

beta_s_comp = [beta_s_true,beta_s_iv, beta_s_gmm]
