%% Import data
%  filename: "C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial
%  organization 1\Part 1 Isis\Homework\data.csv"
%  foldername: 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial
%  organization 1\Part 1 Isis\Homework'
%
% Initialization
clearvars 
% Working Directory
cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Homework'
% Import the data
Data_Initial = csvread("C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Homework\data.csv");


%% Part 1: Estimation of the model

% Key Environmental variable

% Dataset Size
kt = size(Data_Initial,1);



% Variables
year = Data_Initial(:,1);
num_manuf = Data_Initial(:,2);
num_brand = Data_Initial(:,3);
price = Data_Initial(:,4);
X0 = Data_Initial(:,4:7); % Price, Cylinder, Weight, Horsepower
fuel_cost = Data_Initial(:,8);
nb_households = Data_Initial(:,9);
sj = Data_Initial(:,10);
co2 = Data_Initial(:,11);

% Exogenous Variables: Intercept, Price, Cylinder, Weight, Horsepower
X0 = [ones(kt,1) X0];
Z  = [X0 Data_Initial(:,12:23)];
Inverse_ZTZ = inv(Z'*Z);

% Dataset Time Period
T = size(unique(year),1); %T=6

% Computation : Outside Option Market Share

prodsMarket = accumarray(num_market,ones(kt,1));

yearStarts = zeros(T,1);
yearEnds = zeros(T,1);
yearStarts(1) = 1;
yearEnds(1) = prodsYear(1);

for t=2:T
    yearStarts(t) = yearEnds(t-1) + 1;
    yearEnds(t) = yearStarts(t) + prodsYear(t) - 1;
end

s0 = zeros(kt,1);


for tt = 1:T
    Jt = yearStarts(tt):yearEnds(tt);
    s0(Jt) = 1-sum(sj(Jt));
end

%% Part 1: Question 1

% Uniform Draws
rng(2021);
ns = 100;
nu = rand(1,ns);

% Delta initilization

delta_init = log(sj./s0);

% Graphical Exploration of Sigma
% Optimize objective function with sigma over grid [0,1]


grid = [0:0.03:1];
nb_grid = size(grid,2);
obj_value = zeros(nb_grid,1);


for i = 1:nb_grid
    obj_value(i) = gmm_blp(grid(i), sj, X0, fuel_cost, Z, Inverse_ZTZ, nu, delta_init, ...
        yearStarts, yearEnds);
end

% Plot optimized objective function with sigma over grid [0,1]

figure 
plot(grid, obj_value, 'k', 'linewidth', 1.7)
hold on
plot([0.64;0.64], [min(obj_value); max(obj_value)], 'r', 'linewidth', 1.7)% gmm objective function is the lowest around sigma 0.6
xlabel('\sigma')
ylabel('GMM objective function')
axis tight

% Exploration of Sigma with fminunc

opt = optimset('display', 'iter', 'tolX', 1e-4);
sigma_hat_Nfe = fminunc(@(sigma)gmm_blp(sigma, sj, X0, fuel_cost, Z, Inverse_ZTZ, nu, delta_init, ...
    yearStarts, yearEnds), 0.6, opt); % 0.6 initial guess  and sigma=0.6462

% Inversion of the market share using sigma_hat_Nfe
delta_nofe = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_nofe(Jt), pb_cv(tt)] = blp_contraction(sigma_hat_Nfe, sj(Jt), fuel_cost(Jt), delta_init(Jt), nu);   
end

% IV regression without FE
[beta_hat_Nfe] = ivregression(delta_nofe, X0, fuel_cost, Z);

alpha_hat_Nfe = beta_hat_Nfe(size(beta_hat_Nfe,1));

alpha_indi_Nfe =  alpha_hat_Nfe + nu*sigma_hat_Nfe;

%Plot the distribution of the fuel cost sensitivities

[f_alpha_Nfe, x_alpha_Nfe]= ksdensity(alpha_indi_Nfe);

figure
plot(x_alpha_Nfe, f_alpha_Nfe)
xlabel('fuel cost sensitivity')
ylabel('density')


%% Part 1: Question 2

% Construct Brand Dummy

% dummy brand
brand_fe = zeros(kt, max(num_brand));

for t=1:kt
    brand_fe(t,num_brand(t))=1;
end

brand_fe (:, max(num_brand)) = [];

X0_fe = [X0 brand_fe];
Z_fe = [X0_fe Data_Initial(:,12:23)];
Inverse_ZTZ_fe = inv(Z_fe'*Z_fe);


% Graphical Exploration of Sigma
% Optimize objective function with sigma over grid [-2,2]

grid_fe = [-2:0.1:2];
nb_grid_fe = size(grid_fe,2);
obj_value_fe = zeros(nb_grid,1);

for i = 1:nb_grid_fe
    obj_value_fe(i) = gmm_blp(grid_fe(i), sj, X0_fe, fuel_cost, Z_fe, Inverse_ZTZ_fe, nu, delta_init, yearStarts, yearEnds);
end

% Plot optimized objective function with sigma over grid [-2,2]

figure
plot(grid_fe, obj_value_fe) % gmm objective function is the lowest around sigma  -0.5137
hold on
plot([-0.5;-0.5], [min(obj_value_fe); max(obj_value_fe)], 'r',  'linewidth', 1.7)% gmm objective function is the lowest around sigma 0.55
xlabel('\sigma')
ylabel('GMM objective function')
axis tight


% Exploration of Sigma with fminunc
options = optimset('display', 'iter', 'tolX', 1e-4);

sigma_hat_fe = fminunc(@(sigma)gmm_blp(sigma, sj, X0_fe, fuel_cost, Z_fe, Inverse_ZTZ_fe, nu, delta_init, ...
    yearStarts, yearEnds),-0.5, options); % -0.5 initial guess of sigma
 

% Inversion of the market share using sigma_hat_fe
delta_fe = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_fe(Jt), pb_cv_fe(tt)] = blp_contraction(sigma_hat_fe, sj(Jt), fuel_cost(Jt), delta_init(Jt), nu);   
end

% IV regression without FE
[beta_hat_fe] = ivregression(delta_fe, X0_fe, fuel_cost, Z_fe);

alpha_hat_fe = beta_hat_fe(size(beta_hat_fe,1));

alpha_indi_fe =  alpha_hat_fe + nu*sigma_hat_fe;



[f_alpha_fe, x_alpha_fe]= ksdensity(alpha_indi_fe);

%Plot the distribution of the fuel cost sensitivities

figure
plot(x_alpha_fe, f_alpha_fe)
xlabel('fuel cost sensitivity with brand fixed effect')
ylabel('density')

%% Part 1: Question 3


Jt_2008 = yearStarts(T):yearEnds(T);
kt_2008 = size(Jt_2008,2);
sj_2008 = sj(Jt_2008);
num_manuf_2008 = num_manuf(Jt_2008);
price_2008 = price(Jt_2008);
fuel_cost_2008 = fuel_cost(Jt_2008);
delta_fe_2008 = delta_fe(Jt_2008);
delta_fe_2008_np = delta_fe_2008 - beta_hat_fe(2) * price_2008;
co2_2008 = co2(Jt_2008);
% 10 euro for each gram per km 
% remember price unit is in 10,000 euro

tax_2008 = co2(Jt_2008)*10/10000;

numer = exp(delta_fe_2008  + sigma_hat_fe*fuel_cost_2008*nu);
denom = 1+sum(numer,1);
% i for draws
s_hat_ij = numer./denom;
s_hat_j= mean(s_hat_ij,2);  

Temp_omega1 = diag(sum((beta_hat_fe(2)+0).*s_hat_ij,2));
Temp_omega2 = (beta_hat_fe(2)+0).*s_hat_ij*s_hat_ij';

omega = (Temp_omega1-Temp_omega2)/ns;

own = (repmat(num_manuf_2008, 1, kt_2008) == repmat(num_manuf_2008', kt_2008, 1));

omega_tilde = omega.*own;


% marginal cost
mc_2008 = price_2008 + omega_tilde \ sj_2008;
[f_mc_2008, x_mc_2008]= ksdensity(mc_2008);

%Plot the distribution of marginal cost

figure
plot(x_mc_2008, f_mc_2008)
axis([0 max(x_mc_2008) 0 max(1.05*f_mc_2008)])
xlabel('Distribution of Marginal Cost in 2008')
ylabel('density')

Avg_price_2008 = mean(price_2008);
Avg_markup_2008_Def1 = mean((price_2008-mc_2008)./price_2008);
Avg_markup_2008_Def2 = mean((price_2008)./mc_2008);

profit_temp =  (price_2008-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Sum_profit = sum(profit_temp);

Sum_cs = consumer(delta_fe_2008,price_2008 ,beta_hat_fe(2),alpha_hat_fe, fuel_cost(Jt_2008), sigma_hat_fe,nu,ns)*max(nb_households(Jt_2008))/10000;


%% Part 2: Counterfactual Simulations

opt = optimset('Display', 'iter','PlotFcn',@optimplotx ,...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');


price_2008_CF = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, 10 , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, alpha_hat_fe, beta_hat_fe(2), sigma_hat_fe, nu), 1+zeros(kt_2008,1), opt);


% It turns out that new price is the same as old price in 4 digits but
% still different in more digits shown by new == old.
% The influence of such tax is weak 
% when tax rate is 10 euro for each gram per km

price_2008_matrix = [price_2008_CF price_2008 price_2008_CF==price_2008];

Avg_price_2008_CF = mean(price_2008_CF);
Avg_markup_2008_CF_Def1 = mean((price_2008_CF-mc_2008)./price_2008_CF);
Avg_markup_2008_CF_Def2 = mean((price_2008_CF)./mc_2008);

profit_temp_CF =  (price_2008_CF-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Sum_profit_CF = sum(profit_temp_CF);

Sum_cs_CF = consumer(delta_fe_2008_np, price_2008_CF , beta_hat_fe(2) , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe,nu,ns)*max(nb_households(Jt_2008))/10000;

Sum_tax_revenue_CF = tax(price_2008_CF, mc_2008, co2_2008, 10 , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, alpha_hat_fe, beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));
















%%

grid_tax = 0:5:100;
Array_Avg_price_2008_CF = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def1 = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def2 = zeros(1, size(grid_tax,2));
Array_profit_temp_CF = zeros(kt_2008, size(grid_tax,2));
Array_Sum_profit_CF = zeros(1, size(grid_tax,2));
Array_Sum_cs_CF = zeros(1, size(grid_tax,2));
Array_Sum_tax_revenue_CF = zeros(1, size(grid_tax,2));
Matrix_price = zeros(kt_2008, 2*size(grid_tax,2)+1);

opt2 = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

for t = 1:size(grid_tax,2)



price_2008_Temp = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, alpha_hat_fe, beta_hat_fe(2), sigma_hat_fe, nu), 1+zeros(kt_2008,1), opt2);


% It turns out that new price is the same as old price in 4 digits but
% still different in more digits shown by new == old.
% The influence of such tax is weak 
% when tax rate is 10 euro for each gram per km

Matrix_price(:,t) = price_2008_Temp;
Matrix_price(:,t+size(grid_tax,2)) = price_2008_Temp == price_2008;


Array_Avg_price_2008_CF(1,t) = mean(price_2008_Temp);
Array_Avg_markup_2008_CF_Def1(1,t) = mean((price_2008_Temp-mc_2008)./price_2008_Temp);
Array_Avg_markup_2008_CF_Def2(1,t) = mean((price_2008_Temp)./mc_2008);

tax = co2_2008*grid_tax(t)/10000;
numer = exp(delta_fe_2008_np + beta_hat_fe(2)*(price_2008_Temp + tax)  + sigma_hat_fe*fuel_cost_2008*nu);
denom = 1+sum(numer,1);

% i for draws
sij_2008 = numer./denom;

sj_2008= mean(sij_2008,2);  

qj_2008 = max(nb_households(Jt_2008))*sj_2008;

tax_revenue = sum(tax.*qj_2008)*10000/1000000;


Array_profit_temp_CF(:,t) =  (price_2008_Temp-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Array_Sum_profit_CF(1,t) = sum(Array_profit_temp_CF(:,t));

Array_Sum_cs_CF(1,t) = consumer2(delta_fe_2008_np, price_2008_Temp , beta_hat_fe(2) , tax , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe, nu , ns)*max(nb_households(Jt_2008))/10000;

% Array_Sum_tax_revenue_CF(1,t) = tax(price_2008_Temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
%    beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));

    Array_Sum_tax_revenue_CF(1,t) = tax_revenue;
end

Matrix_price(:,2*size(grid_tax,2)+1) = price_2008;


figure
plot(grid_tax,Array_Sum_tax_revenue_CF);
hold on
plot(grid_tax,Array_Sum_profit_CF);
plot(grid_tax,Array_Sum_cs_CF);
legend('Tax', 'Profit', 'ConsumerSurplus')
axis tight
hold off


%% A quick look into results of simulation without waiting for loops
% 
% load test2.mat;
% 
% for t = 1:size(grid_tax,2)
% 
% 
% 
% price_2008_Temp = Matrix_price(:,t);
% 
% % It turns out that new price is the same as old price in 4 digits but
% % still different in more digits shown by new == old.
% % The influence of such tax is weak 
% % when tax rate is 10 euro for each gram per km
% 
% 
% Array_Avg_price_2008_CF(1,t) = mean(price_2008_Temp);
% Array_Avg_markup_2008_CF_Def1(1,t) = mean((price_2008_Temp-mc_2008)./price_2008_Temp);
% Array_Avg_markup_2008_CF_Def2(1,t) = mean((price_2008_Temp)./mc_2008);
% 
% tax = co2_2008*grid_tax(t)/10000;
% numer = exp(delta_fe_2008_np + beta_hat_fe(2)*(price_2008_Temp + tax)  + sigma_hat_fe*fuel_cost_2008*nu);
% denom = 1+sum(numer,1);
% 
% % i for draws
% sij_2008 = numer./denom;
% 
% sj_2008= mean(sij_2008,2);  
% 
% qj_2008 = max(nb_households(Jt_2008))*sj_2008;
% 
% tax_revenue = sum(tax.*qj_2008)*10000/1000000;
% 
% 
% Array_profit_temp_CF(:,t) =  (price_2008_Temp-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 
% 
% Array_Sum_profit_CF(1,t) = sum(Array_profit_temp_CF(:,t));
% 
% Array_Sum_cs_CF(1,t) = consumer2(delta_fe_2008_np, price_2008_Temp , beta_hat_fe(2) , tax , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe, nu , ns)*max(nb_households(Jt_2008))/10000;
% 
% % Array_Sum_tax_revenue_CF(1,t) = tax(price_2008_Temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
% %    beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));
% 
%     Array_Sum_tax_revenue_CF(1,t) = tax_revenue;
% end
% 
% figure
% plot(grid_tax,Array_Sum_tax_revenue_CF);
% hold on
% plot(grid_tax,Array_Sum_profit_CF);
% plot(grid_tax,Array_Sum_cs_CF);
% legend('Tax', 'Profit', 'ConsumerSurplus')
% axis tight
% hold off


%% test on grid [0,10000]

grid_tax = 0:1000:10000;
Array_Avg_price_2008_CF = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def1 = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def2 = zeros(1, size(grid_tax,2));
Array_profit_temp_CF = zeros(kt_2008, size(grid_tax,2));
Array_Sum_profit_CF = zeros(1, size(grid_tax,2));
Array_Sum_cs_CF = zeros(1, size(grid_tax,2));
Array_Sum_tax_revenue_CF = zeros(1, size(grid_tax,2));
Matrix_price = zeros(kt_2008, 2*size(grid_tax,2)+1);

opt2 = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

for t = 1:size(grid_tax,2)



price_2008_Temp = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, alpha_hat_fe, beta_hat_fe(2), sigma_hat_fe, nu), price_2008, opt2);


% It turns out that new price is the same as old price in 4 digits but
% still different in more digits shown by new == old.
% The influence of such tax is weak 
% when tax rate is 10 euro for each gram per km

Matrix_price(:,t) = price_2008_Temp;
Matrix_price(:,t+size(grid_tax,2)) = price_2008_Temp == price_2008;


Array_Avg_price_2008_CF(1,t) = mean(price_2008_Temp);
Array_Avg_markup_2008_CF_Def1(1,t) = mean((price_2008_Temp-mc_2008)./price_2008_Temp);
Array_Avg_markup_2008_CF_Def2(1,t) = mean((price_2008_Temp)./mc_2008);

tax = co2_2008*grid_tax(t)/10000;
numer = exp(delta_fe_2008_np + beta_hat_fe(2)*(price_2008_Temp + tax)  + sigma_hat_fe*fuel_cost_2008*nu);
denom = 1+sum(numer,1);

% i for draws
sij_2008 = numer./denom;

sj_2008= mean(sij_2008,2);  

qj_2008 = max(nb_households(Jt_2008))*sj_2008;

tax_revenue = sum(tax.*qj_2008)*10000/1000000;


Array_profit_temp_CF(:,t) =  (price_2008_Temp-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Array_Sum_profit_CF(1,t) = sum(Array_profit_temp_CF(:,t));

Array_Sum_cs_CF(1,t) = consumer2(delta_fe_2008_np, price_2008_Temp , beta_hat_fe(2) , tax , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe, nu , ns)*max(nb_households(Jt_2008))/10000;

% Array_Sum_tax_revenue_CF(1,t) = tax(price_2008_Temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
%    beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));

    Array_Sum_tax_revenue_CF(1,t) = tax_revenue;
end

Matrix_price(:,2*size(grid_tax,2)+1) = price_2008;


figure
plot(grid_tax,Array_Sum_tax_revenue_CF);
hold on
plot(grid_tax,Array_Sum_profit_CF);
plot(grid_tax,Array_Sum_cs_CF);

save 0-10000test.mat
legend('Tax', 'Profit', 'ConsumerSurplus')
axis tight
hold off

%% test on grid [0,2000]

grid_tax = 0:200:2000;
Array_Avg_price_2008_CF = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def1 = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def2 = zeros(1, size(grid_tax,2));
Array_profit_temp_CF = zeros(kt_2008, size(grid_tax,2));
Array_Sum_profit_CF = zeros(1, size(grid_tax,2));
Array_Sum_cs_CF = zeros(1, size(grid_tax,2));
Array_Sum_tax_revenue_CF = zeros(1, size(grid_tax,2));
Matrix_price = zeros(kt_2008, 2*size(grid_tax,2)+1);

opt2 = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

for t = 1:size(grid_tax,2)



price_2008_Temp = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, alpha_hat_fe, beta_hat_fe(2), sigma_hat_fe, nu), price_2008, opt2);


% It turns out that new price is the same as old price in 4 digits but
% still different in more digits shown by new == old.
% The influence of such tax is weak 
% when tax rate is 10 euro for each gram per km

Matrix_price(:,t) = price_2008_Temp;
Matrix_price(:,t+size(grid_tax,2)) = price_2008_Temp == price_2008;


Array_Avg_price_2008_CF(1,t) = mean(price_2008_Temp);
Array_Avg_markup_2008_CF_Def1(1,t) = mean((price_2008_Temp-mc_2008)./price_2008_Temp);
Array_Avg_markup_2008_CF_Def2(1,t) = mean((price_2008_Temp)./mc_2008);

tax = co2_2008*grid_tax(t)/10000;
numer = exp(delta_fe_2008_np + beta_hat_fe(2)*(price_2008_Temp + tax)  + sigma_hat_fe*fuel_cost_2008*nu);
denom = 1+sum(numer,1);

% i for draws
sij_2008 = numer./denom;

sj_2008= mean(sij_2008,2);  

qj_2008 = max(nb_households(Jt_2008))*sj_2008;

tax_revenue = sum(tax.*qj_2008)*10000/1000000;


Array_profit_temp_CF(:,t) =  (price_2008_Temp-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Array_Sum_profit_CF(1,t) = sum(Array_profit_temp_CF(:,t));

Array_Sum_cs_CF(1,t) = consumer2(delta_fe_2008_np, price_2008_Temp , beta_hat_fe(2) , tax , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe, nu , ns)*max(nb_households(Jt_2008))/10000;

% Array_Sum_tax_revenue_CF(1,t) = tax(price_2008_Temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
%    beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));

    Array_Sum_tax_revenue_CF(1,t) = tax_revenue;
end

Matrix_price(:,2*size(grid_tax,2)+1) = price_2008;


figure
plot(grid_tax,Array_Sum_tax_revenue_CF);
hold on
plot(grid_tax,Array_Sum_profit_CF);
plot(grid_tax,Array_Sum_cs_CF);
legend('Tax', 'Profit', 'ConsumerSurplus')
axis tight
hold off

save 0-2000test.mat
legend('Tax', 'Profit', 'ConsumerSurplus')
axis tight
hold off


%% test on grid [0,1400]

grid_tax = 0:20:1400;
Array_Avg_price_2008_CF = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def1 = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def2 = zeros(1, size(grid_tax,2));
Array_profit_temp_CF = zeros(kt_2008, size(grid_tax,2));
Array_Sum_profit_CF = zeros(1, size(grid_tax,2));
Array_Sum_cs_CF = zeros(1, size(grid_tax,2));
Array_Sum_tax_revenue_CF = zeros(1, size(grid_tax,2));
Matrix_price = zeros(kt_2008, 2*size(grid_tax,2)+1);

opt2 = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

for t = 1:size(grid_tax,2)



price_2008_Temp = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, alpha_hat_fe, beta_hat_fe(2), sigma_hat_fe, nu), price_2008, opt2);


% It turns out that new price is the same as old price in 4 digits but
% still different in more digits shown by new == old.
% The influence of such tax is weak 
% when tax rate is 10 euro for each gram per km

Matrix_price(:,t) = price_2008_Temp;
Matrix_price(:,t+size(grid_tax,2)) = price_2008_Temp == price_2008;


Array_Avg_price_2008_CF(1,t) = mean(price_2008_Temp);
Array_Avg_markup_2008_CF_Def1(1,t) = mean((price_2008_Temp-mc_2008)./price_2008_Temp);
Array_Avg_markup_2008_CF_Def2(1,t) = mean((price_2008_Temp)./mc_2008);

tax = co2_2008*grid_tax(t)/10000;
numer = exp(delta_fe_2008_np + beta_hat_fe(2)*(price_2008_Temp + tax)  + sigma_hat_fe*fuel_cost_2008*nu);
denom = 1+sum(numer,1);

% i for draws
sij_2008 = numer./denom;

sj_2008= mean(sij_2008,2);  

qj_2008 = max(nb_households(Jt_2008))*sj_2008;

tax_revenue = sum(tax.*qj_2008)*10000/1000000;


Array_profit_temp_CF(:,t) =  (price_2008_Temp-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Array_Sum_profit_CF(1,t) = sum(Array_profit_temp_CF(:,t));

Array_Sum_cs_CF(1,t) = consumer2(delta_fe_2008_np, price_2008_Temp , beta_hat_fe(2) , tax , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe, nu , ns)*max(nb_households(Jt_2008))/10000;

% Array_Sum_tax_revenue_CF(1,t) = tax(price_2008_Temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
%    beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));

    Array_Sum_tax_revenue_CF(1,t) = tax_revenue;
end

Matrix_price(:,2*size(grid_tax,2)+1) = price_2008;


figure
plot(grid_tax,Array_Sum_tax_revenue_CF);
hold on
plot(grid_tax,Array_Sum_profit_CF);
plot(grid_tax,Array_Sum_cs_CF);
legend('Tax', 'Profit', 'ConsumerSurplus')
axis tight
hold off

save 0-1400test.mat
legend('Tax', 'Profit', 'ConsumerSurplus')
axis tight
hold off
