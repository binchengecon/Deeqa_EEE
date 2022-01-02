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
Data_Initial = csvread("data.csv");


%% Part 1: Estimation of the model

% Key Environmental variable

% Dataset Size
kt = size(Data_Initial,1);



% Variables
year = Data_Initial(:,1);
num_manuf = Data_Initial(:,2);
num_brand = Data_Initial(:,3);
price = Data_Initial(:,4);

X0 = Data_Initial(:,5:7); % ,Cylinder, Weight, Horsepower

fuel_cost = Data_Initial(:,8);
nb_households = Data_Initial(:,9);
sj = Data_Initial(:,10);
co2 = Data_Initial(:,11);

% Exogenous Variables: Intercept, Cylinder, Weight, Horsepower
X0 = [ones(kt,1) X0 fuel_cost];

% Endogenous Variables: Price

Endo = price;


Z  = [X0 Data_Initial(:,12:23)];
Inverse_ZTZ = inv(Z'*Z);

% Dataset Time Period
T = size(unique(year),1); %T=6

% Computation : Outside Option Market Share
num_market = year- 2002;
prodsMarket = accumarray(num_market,ones(kt,1));

yearStarts = zeros(T,1);
yearEnds = zeros(T,1);
yearStarts(1) = 1;
yearEnds(1) = prodsMarket(1);

for t=2:T
    yearStarts(t) = yearEnds(t-1) + 1;
    yearEnds(t) = yearStarts(t) + prodsMarket(t) - 1;
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
nu = rand(1,ns)-1/2; 
%U[-1/2,1/2]
% nu = rand(1,ns);
%U[0,1]

% Delta initilization

delta_init = log(sj./s0);

% Graphical Exploration of Sigma
% Optimize objective function with sigma over grid [-4,4]

grid = [-4:0.1:4];

% grid = [0:0.1:4];
nb_grid = size(grid,2);
obj_value = zeros(nb_grid,1);

% grid(i) represent sigma
% obj_value represent the gmm value under sigma : gmm_blp(grid(i),...)

for i = 1:nb_grid
    obj_value(i) = gmm_blp(grid(i), sj, X0, fuel_cost, Endo, Z, Inverse_ZTZ, nu, delta_init, yearStarts, yearEnds);
end

% % Plot optimized objective function with sigma over grid [-4,4]
% 
[~,Index_min_obj] = min(obj_value);
% 
% figure 
% plot(grid, obj_value)
% hold on
% plot([grid(Index_min_obj);grid(Index_min_obj)], [min(obj_value); max(obj_value)], 'r', 'linewidth', 1.7)
% % gmm objective function is the lowest around sigma -0.3122
% xlabel('\sigma')
% ylabel('GMM objective function')
% title("\sigma =" + sigma_hat_Nfe)
% axis tight
% 
% saveas(gcf, 'GMM_NfeNegativeSigma.jpg' );

% Exploration of Sigma with fminunc

opt = optimset('display', 'iter', 'tolX', 1e-4);
sigma_hat_Nfe = fminunc(@(sigma)gmm_blp(sigma, sj, X0, fuel_cost, Endo, Z, Inverse_ZTZ, nu, delta_init, ...
    yearStarts, yearEnds), grid(Index_min_obj), opt); % grid(Index_min_obj) initial guess  and sigma= sigma_hat_Nfe

% Inversion of the market share using sigma_hat_Nfe

delta_nofe = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_nofe(Jt), pb_cv(tt)] = blp_contraction(sigma_hat_Nfe, sj(Jt), delta_init(Jt), nu,fuel_cost(Jt,:));   
end


% IV regression without FE
[beta_hat_Nfe] = ivregression(delta_nofe, X0, Endo, Z);

% fuel cost is the last varibale in exogenous varibale

alpha_hat_Nfe = beta_hat_Nfe(end-1);

alpha_indi_Nfe =  alpha_hat_Nfe + nu*sigma_hat_Nfe;

% %Plot the distribution of the fuel cost sensitivities
% 
% [f_alpha_Nfe, x_alpha_Nfe]= ksdensity(alpha_indi_Nfe);
% 
% figure
% plot(x_alpha_Nfe, f_alpha_Nfe)
% xlabel('fuel cost sensitivity')
% ylabel('density')


%% Part 1: Question 2

% Construct Brand Dummy
% dummy year

year_fe = zeros(kt,T);
for tt = 1:T
    Jt = yearStarts(tt):yearEnds(tt);
    year_fe(Jt,tt)=1;
end
year_fe (:, 1) = [];

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
% Optimize objective function with sigma over grid [-4,4]

grid_fe = [-4:0.1:4];

% grid_fe = [0:0.1:4];
nb_grid_fe = size(grid_fe,2);
obj_value_fe = zeros(nb_grid_fe,1);

for i = 1:nb_grid_fe
    obj_value_fe(i) = gmm_blp(grid_fe(i), sj, X0_fe, fuel_cost, Endo, Z_fe, Inverse_ZTZ_fe, nu, delta_init, yearStarts, yearEnds);
end

% % Plot optimized objective function with sigma over grid [-4,4]
% 
[~,Index_min_obj_fe] = min(obj_value_fe);
% 
% figure
% plot(grid_fe, obj_value_fe) % gmm objective function is the lowest around sigma  -0.4030
% hold on
% plot([grid_fe(Index_min_obj_fe);grid_fe(Index_min_obj_fe)], [min(obj_value_fe); max(obj_value_fe)], 'r',  'linewidth', 1.7)% gmm objective function is the lowest around sigma -0.4030
% xlabel('\sigma')
% ylabel('GMM objective function')
% title("\sigma =" + sigma_hat_fe)
% axis tight
% 
% saveas(gcf, 'GMM_feNegativeSigma.jpg' );

% Exploration of Sigma with fminunc
options = optimset('display', 'iter', 'tolX', 1e-4);

sigma_hat_fe = fminunc(@(sigma)gmm_blp(sigma, sj, X0_fe, fuel_cost, Endo, Z_fe, Inverse_ZTZ_fe, nu, delta_init, ...
    yearStarts, yearEnds),grid_fe(Index_min_obj_fe), options); % -0.4 initial guess, sigma = -0.4030
 

% Inversion of the market share using sigma_hat_fe
delta_fe = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_fe(Jt), pb_cv_fe(tt)] = blp_contraction(sigma_hat_fe, sj(Jt), delta_init(Jt), nu,fuel_cost(Jt,:));   
end

% IV regression without FE
[beta_hat_fe] = ivregression(delta_fe, X0_fe, Endo, Z_fe);

% fuel_cost, brand
alpha_hat_fe = beta_hat_fe(5);

alpha_indi_fe =  alpha_hat_fe + nu*sigma_hat_fe;

% 
% 
% 
% %Plot the distribution of the fuel cost sensitivities
% [f_alpha_fe, x_alpha_fe]= ksdensity(alpha_indi_fe);
% 
% figure
% plot(x_alpha_fe, f_alpha_fe)
% xlabel('fuel cost sensitivity with brand fixed effect')
% ylabel('density')
% 
% 
% figure 
% plot(x_alpha_fe, f_alpha_fe)
% hold on
% plot(x_alpha_Nfe, f_alpha_Nfe)
% xlabel('Fuel cost sensitivity')
% ylabel('Density')
% legend('With brand fe','Without brand fe')
% title("\sigma =" + sigma_hat_fe)
% hold off
% 
% saveas(gcf,'sensitivityNegativeSigma.jpg');


%% Part 1: Question 2.5


X0_fe2 = [X0 brand_fe year_fe];
Z_fe2 = [X0_fe2 Data_Initial(:,12:23)];
Inverse_ZTZ_fe2 = inv(Z_fe2'*Z_fe2);


% Graphical Exploration of Sigma
% Optimize objective function with sigma over grid [-4,4]

grid_fe2 = [-10:0.5:10];

% grid_fe = [0:0.1:4];
nb_grid_fe2 = size(grid_fe2,2);
obj_value_fe2 = zeros(nb_grid_fe2,1);

for i = 1:nb_grid_fe2
    obj_value_fe2(i) = gmm_blp(grid_fe2(i), sj, X0_fe2, fuel_cost, Endo, Z_fe2, Inverse_ZTZ_fe2, nu, delta_init, yearStarts, yearEnds);
end

% % Plot optimized objective function with sigma over grid [-4,4]
% 
[~,Index_min_obj_fe2] = min(obj_value_fe2);
% 
% figure
% plot(grid_fe, obj_value_fe) % gmm objective function is the lowest around sigma  -0.4030
% hold on
% plot([grid_fe(Index_min_obj_fe);grid_fe(Index_min_obj_fe)], [min(obj_value_fe); max(obj_value_fe)], 'r',  'linewidth', 1.7)% gmm objective function is the lowest around sigma -0.4030
% xlabel('\sigma')
% ylabel('GMM objective function')
% title("\sigma =" + sigma_hat_fe)
% axis tight
% 
% saveas(gcf, 'GMM_feNegativeSigma.jpg' );

% Exploration of Sigma with fminunc
options = optimset('display', 'iter', 'tolX', 1e-4);

sigma_hat_fe2 = fminunc(@(sigma)gmm_blp(sigma, sj, X0_fe2, fuel_cost, Endo, Z_fe2, Inverse_ZTZ_fe2, nu, delta_init, ...
    yearStarts, yearEnds),grid_fe2(Index_min_obj_fe2), options); % -0.4 initial guess, sigma = -0.4030
 

% Inversion of the market share using sigma_hat_fe
delta_fe2 = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_fe2(Jt), pb_cv_fe2(tt)] = blp_contraction(sigma_hat_fe2, sj(Jt), delta_init(Jt), nu,fuel_cost(Jt,:));   
end

% IV regression without FE
[beta_hat_fe2] = ivregression(delta_fe2, X0_fe2, Endo, Z_fe2);

% fuel_cost, brand
alpha_hat_fe2 = beta_hat_fe2(5);

alpha_indi_fe2 =  alpha_hat_fe2 + nu*sigma_hat_fe2;

%% Part 1: Question 2.8

% Collect data into tables for exportation into Latex


Table_Regression = zeros(7,3);

% No FE case
% beta_hat_Nfe: 1, cy,w,hp,    fuelcost, price

Table_Regression(1:5,1) = beta_hat_Nfe(1:5);
Table_Regression(6,1)   = sigma_hat_Nfe;
Table_Regression(7,1)   = beta_hat_Nfe(end);

%FE case with brand effect
% beta_hat_fe: 1, cy,w,hp, fuelcost,   brand(40 brands, 39 brand_fe) , price
%              1  2-4         5          6-44,                          45

Table_Regression(1:5,2) = beta_hat_fe(1:5);
Table_Regression(6,2) = sigma_hat_fe;
Table_Regression(7,2)   = beta_hat_fe(end);

%FE case with brand and year effect
% beta_hat_fe2: 1, cy,w,hp, fuelcost,   brand(40 brands, 39 brand_fe) , year (6 years, 5 year_fe), price
%              1    2-4       5          6-44,                          45-49                      50

Table_Regression(1:5,3) = beta_hat_fe2(1:5);
Table_Regression(6,3) = sigma_hat_fe2;
Table_Regression(7,3)   = beta_hat_fe2(end);



ReResultVarNames = {'No Fixed Effect','With Brand Fixed Effect','With Brand and Year Fixed Effect'};
% Table_CFResult = renamevars(Table_CFResult,allvars,CFResultNames);
ReResultRowNames = {'Intercept','Cylinder','Weight','Horsepower','Fuel Cost','Sigma','Price'};

Table_Regression   = array2table(Table_Regression,'VariableNames',ReResultVarNames,'RowNames',ReResultRowNames);


%% Part 1: Question 3


Jt_2008 = yearStarts(T):yearEnds(T);
kt_2008 = size(Jt_2008,2);
sj_2008 = sj(Jt_2008);
num_manuf_2008 = num_manuf(Jt_2008);
price_2008 = price(Jt_2008);
fuel_cost_2008 = fuel_cost(Jt_2008);
Endo_2008 = Endo(Jt_2008,:);
delta_fe_2008 = delta_fe(Jt_2008);


% Intercept, Exo_Except Fuel cost, fuel cost, brand fixed effect, price
%     1          2-4                    5                          end

delta_fe_2008_np = delta_fe_2008 - beta_hat_fe(end) * price_2008;
co2_2008 = co2(Jt_2008);

% 10 euro for each gram per km 
% remember price unit is in 10,000 euro
% 140 g/km -> 140*10 euros->1400/10000 in 10,000 euros
% comparable with price

% Temp_omega1 = diag(sum((beta_hat_fe(end-1)+zeros(1,ns)).*s_hat_ij,2));
% Temp_omega2 = (beta_hat_fe(end-1)+zeros(1,ns)).*s_hat_ij*s_hat_ij';

% omega = (1/ns)*(Temp_omega1-Temp_omega2);

% own = (repmat(num_manuf_2008, 1, kt_2008) == repmat(num_manuf_2008', kt_2008, 1));
% 
% omega_tilde = omega.*own;

tax_2008 = co2_2008*10/10000;

numer = exp(delta_fe_2008  + sigma_hat_fe*fuel_cost_2008*nu);

denom = 1+sum(numer,1);

% i for draws
s_hat_ij = numer./denom;
s_hat_j= mean(s_hat_ij,2);  



Temp_omega1 = diag(sum((beta_hat_fe(end)).*s_hat_j,2));
Temp_omega2 = (beta_hat_fe(end)).*s_hat_j*s_hat_j';

omega = (Temp_omega1-Temp_omega2);

own = (repmat(num_manuf_2008, 1, kt_2008) == repmat(num_manuf_2008', kt_2008, 1));

omega_tilde = omega.*own;


% marginal cost
mc_2008 = price_2008 + omega_tilde \ sj_2008;

% %Plot the distribution of marginal cost
% [f_mc_2008, x_mc_2008]= ksdensity(mc_2008,'Support','positive','BoundaryCorrection','reflection');

% figure
% plot(x_mc_2008, f_mc_2008)
% axis([min(x_mc_2008) max(x_mc_2008) 0 max(1.05*f_mc_2008)])
% xlabel('Distribution of Marginal Cost in 2008')
% ylabel('density')
% title("\sigma =" + sigma_hat_fe)
% 
% 
% saveas(gcf,'MarginalCost2008NegativeSigma.jpg');

Avg_price_2008 = mean(price_2008);
Avg_markup_2008_Def1 = mean((price_2008-mc_2008)./price_2008);
Avg_markup_2008_Def2 = mean((price_2008)./mc_2008-1);

profit_temp =  (price_2008-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Sum_profit = sum(profit_temp);

Sum_cs = consumer(delta_fe_2008,alpha_hat_fe, fuel_cost_2008, sigma_hat_fe,nu,ns)*max(nb_households(Jt_2008))*10000/1000000;



%% Part 2: Counterfactual Simulations

opt = optimset('Display', 'iter','PlotFcn',@optimplotx ,...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

opt2 = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

opt3 = optimset('Display', 'none',...
    'TolFun',1e-4,'TolX',1e-6,...
    'MaxFunEvals', 1000000,'MaxIter',5000,...
    'LargeScale','off');

price_2008_CF = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, 10 , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, beta_hat_fe(end), sigma_hat_fe, nu), price_2008, opt2);


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

Sum_cs_CF = consumer2(delta_fe_2008_np, price_2008_CF , beta_hat_fe(end) , tax_2008, alpha_hat_fe,fuel_cost_2008,sigma_hat_fe,nu,ns)*max(nb_households(Jt_2008))/10000;

Sum_tax_revenue_CF = tax(price_2008_CF, co2_2008, 10 , fuel_cost_2008, delta_fe_2008_np, ...
      beta_hat_fe(end), sigma_hat_fe, nu,max(nb_households(Jt_2008)));


%% Part 2: Question 2

grid_tax = 0:1:100;
Array_Avg_price_2008_CF = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def1 = zeros(1, size(grid_tax,2));
Array_Avg_markup_2008_CF_Def2 = zeros(1, size(grid_tax,2));
Array_profit_temp_CF = zeros(kt_2008, size(grid_tax,2));
Array_Sum_profit_CF = zeros(1, size(grid_tax,2));
Array_Sum_cs_CF = zeros(1, size(grid_tax,2));
Array_Sum_tax_revenue_CF = zeros(1, size(grid_tax,2));
Matrix_price = zeros(kt_2008, 2*size(grid_tax,2)+1);



for t = 1:size(grid_tax,2)

Progress = grid_tax(t)/max(grid_tax)*100;
Message = sprintf('The loop is running at %0.2f %%',Progress);

price_2008_Temp = fsolve(@(price_temp) CounterFactualSystem(price_temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
    num_manuf_2008, beta_hat_fe(end), sigma_hat_fe, nu), price_2008, opt3);


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
numer = exp(delta_fe_2008_np + beta_hat_fe(end)*(price_2008_Temp + tax)  + sigma_hat_fe*fuel_cost_2008*nu);
denom = 1+sum(numer,1);

% i for draws
sij_2008 = numer./denom;

sj_2008= mean(sij_2008,2);  

qj_2008 = max(nb_households(Jt_2008))*sj_2008;

tax_revenue = sum(tax.*qj_2008)*10000/1000000;


Array_profit_temp_CF(:,t) =  (price_2008_Temp-mc_2008).*sj_2008.*nb_households(Jt_2008)*10000/1000000; 

Array_Sum_profit_CF(1,t) = sum(Array_profit_temp_CF(:,t));

Array_Sum_cs_CF(1,t) = consumer2(delta_fe_2008_np, price_2008_Temp , beta_hat_fe(end) , tax , alpha_hat_fe,fuel_cost_2008,sigma_hat_fe, nu , ns)*max(nb_households(Jt_2008))/10000;

% Array_Sum_tax_revenue_CF(1,t) = tax(price_2008_Temp, mc_2008, co2_2008, grid_tax(t) , fuel_cost_2008, delta_fe_2008_np, ...
%    beta_hat_fe(2), sigma_hat_fe, nu,max(nb_households(Jt_2008)));

    Array_Sum_tax_revenue_CF(1,t) = tax_revenue;
    
    display(Message)

end

Matrix_price(:,2*size(grid_tax,2)+1) = price_2008;


% figure
% plot(grid_tax,Array_Sum_tax_revenue_CF);
% hold on
% plot(grid_tax,Array_Sum_profit_CF);
% plot(grid_tax,Array_Sum_cs_CF);
% legend('Tax Revenue', 'Profit', 'ConsumerSurplus')
% xlabel('Tax Rate')
% ylabel('Amount')
% title("\sigma =" + sigma_hat_fe)
% axis tight
% hold off
% 
% saveas(gcf,'Tax_RevenueNegativeSigma.jpg')

% save changeendNegativeSigma.mat

%% Part 2: Question 2.5

Table_CFResult   = zeros(7,4);
% t=0
Table_CFResult(1,1) = mean(mc_2008);
Table_CFResult(2,1) = Array_Avg_price_2008_CF(1);
Table_CFResult(3,1) = Array_Avg_markup_2008_CF_Def1(1);
Table_CFResult(4,1) = Array_Avg_markup_2008_CF_Def2(1)-1;
Table_CFResult(5,1) = Array_Sum_profit_CF(1);
Table_CFResult(6,1) = Array_Sum_cs_CF(1);
Table_CFResult(7,1) = Array_Sum_tax_revenue_CF(1);


% t = 10
Table_CFResult(1,2) = mean(mc_2008);
Table_CFResult(2,2) = Array_Avg_price_2008_CF(11);
Table_CFResult(3,2) = Array_Avg_markup_2008_CF_Def1(11);
Table_CFResult(4,2) = Array_Avg_markup_2008_CF_Def2(11)-1;
Table_CFResult(5,2) = Array_Sum_profit_CF(11);
Table_CFResult(6,2) = Array_Sum_cs_CF(11);
Table_CFResult(7,2) = Array_Sum_tax_revenue_CF(11);

% t=
[~,Index_Max_TaxR] = max(Array_Sum_tax_revenue_CF);
taxrate_optimal    = grid_tax(Index_Max_TaxR);

Table_CFResult(1,3) = mean(mc_2008);
Table_CFResult(2,3) = Array_Avg_price_2008_CF(Index_Max_TaxR);
Table_CFResult(3,3) = Array_Avg_markup_2008_CF_Def1(Index_Max_TaxR);
Table_CFResult(4,3) = Array_Avg_markup_2008_CF_Def2(Index_Max_TaxR)-1;
Table_CFResult(5,3) = Array_Sum_profit_CF(Index_Max_TaxR);
Table_CFResult(6,3) = Array_Sum_cs_CF(Index_Max_TaxR);
Table_CFResult(7,3) = Array_Sum_tax_revenue_CF(Index_Max_TaxR);

% t=100

Table_CFResult(1,4) = mean(mc_2008);
Table_CFResult(2,4) = Array_Avg_price_2008_CF(100);
Table_CFResult(3,4) = Array_Avg_markup_2008_CF_Def1(100);
Table_CFResult(4,4) = Array_Avg_markup_2008_CF_Def2(100)-1;
Table_CFResult(5,4) = Array_Sum_profit_CF(100);
Table_CFResult(6,4) = Array_Sum_cs_CF(100);
Table_CFResult(7,4) = Array_Sum_tax_revenue_CF(100);


Tax_rate = [0,10,taxrate_optimal,100];
allvars = 1:width(Table_CFResult);
CFResultVarNames = append("Tax Rate =",string(Tax_rate));
% Table_CFResult = renamevars(Table_CFResult,allvars,CFResultNames);
CFResultRowNames = {'Marginal Cost','Average Price', 'Markup (1-c/p)','Markup (p/c-1)','Total Profit', 'Total Consumer Surplus', 'Total Tax Revenue'};

Table_CFResult   = array2table(Table_CFResult,'VariableNames',CFResultVarNames,'RowNames',CFResultRowNames);

%% A quick look into results of simulation without waiting for loops
% 
load changeendNegativeSigma.mat;
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

%% Figure

% Plot optimized objective function with sigma_Nfe over grid [-4,4]

[~,Index_min_obj] = min(obj_value);

figure 
plot(grid, obj_value)
hold on
plot([grid(Index_min_obj);grid(Index_min_obj)], [min(obj_value); max(obj_value)], 'r', 'linewidth', 1.7)
% gmm objective function is the lowest around sigma -0.3122
xlabel('\sigma')
ylabel('GMM objective function')
title("\sigma =" + sigma_hat_Nfe)
axis tight

saveas(gcf, 'GMM_NfeNegativeSigma.jpg' );
saveas(gcf, 'GMM_NfeNegativeSigma.pdf' );

%Plot the distribution of the fuel cost sensitivities

[f_alpha_Nfe, x_alpha_Nfe]= ksdensity(alpha_indi_Nfe);

% figure
% plot(x_alpha_Nfe, f_alpha_Nfe)
% xlabel('fuel cost sensitivity')
% ylabel('density')

% Plot optimized objective function with sigma_fe over grid [-4,4]

[~,Index_min_obj_fe] = min(obj_value_fe);

figure
plot(grid_fe, obj_value_fe) % gmm objective function is the lowest around sigma  -0.4030
hold on
plot([grid_fe(Index_min_obj_fe);grid_fe(Index_min_obj_fe)], [min(obj_value_fe); max(obj_value_fe)], 'r',  'linewidth', 1.7)% gmm objective function is the lowest around sigma -0.4030
xlabel('\sigma')
ylabel('GMM objective function')
title("\sigma =" + sigma_hat_fe)
axis tight

saveas(gcf, 'GMM_feNegativeSigma.jpg' );
saveas(gcf, 'GMM_feNegativeSigma.pdf' );



% Plot the distribution of the fuel cost sensitivities_fe
[f_alpha_fe, x_alpha_fe]= ksdensity(alpha_indi_fe);
[f_alpha_fe2, x_alpha_fe2]= ksdensity(alpha_indi_fe2);

% figure
% plot(x_alpha_fe, f_alpha_fe)
% xlabel('fuel cost sensitivity with brand fixed effect')
% ylabel('density')

% Plot the comparison of fuel cost sensitivity between Nfe and fe


figure 
plot(x_alpha_fe, f_alpha_fe)
hold on
plot(x_alpha_Nfe, f_alpha_Nfe)
xlabel('Fuel cost sensitivity')
ylabel('Density')
legend('With brand fe','Without brand fe')
title("\sigma =" + sigma_hat_fe)
hold off

saveas(gcf,'sensitivityNegativeSigma.jpg');
saveas(gcf,'sensitivityNegativeSigma.pdf');

figure 
plot(x_alpha_fe, f_alpha_fe)
hold on
plot(x_alpha_Nfe, f_alpha_Nfe)
plot(x_alpha_fe2, f_alpha_fe2)
xlabel('Fuel cost sensitivity')
ylabel('Density')
legend('With brand fe','Without brand fe','with brand and year fe')
title("\sigma_b =" + sigma_hat_fe +"    \sigma_{by} =" + sigma_hat_fe2)
hold off

saveas(gcf,'sensitivityNegativeSigma2.jpg');
saveas(gcf,'sensitivityNegativeSigma2.pdf');

%Plot the distribution of marginal cost

[f_mc_2008, x_mc_2008]= ksdensity(mc_2008,'Support','positive','BoundaryCorrection','reflection');

figure
plot(x_mc_2008, f_mc_2008)
axis([min(x_mc_2008) max(x_mc_2008) 0 max(1.05*f_mc_2008)])
xlabel('Distribution of Marginal Cost in 2008')
ylabel('density')
title("\sigma =" + sigma_hat_fe)


saveas(gcf,'MarginalCost2008NegativeSigma.jpg');
saveas(gcf,'MarginalCost2008NegativeSigma.pdf');

%Plot the distribution of tax,profit and consumer surplus

figure
plot(grid_tax,Array_Sum_tax_revenue_CF);
hold on
plot(grid_tax,Array_Sum_profit_CF);
plot(grid_tax,Array_Sum_cs_CF);
legend('Tax Revenue', 'Profit', 'ConsumerSurplus')
xlabel('Tax Rate')
ylabel('Amount')
title("\sigma =" + sigma_hat_fe)
axis tight
hold off

saveas(gcf,'Tax_RevenueNegativeSigma.jpg')
saveas(gcf,'Tax_RevenueNegativeSigma.pdf')

table2latex(Table_CFResult,'Table_CFResult.tex');
table2latex(Table_Regression,'Table_Regression.tex');

