%% Initialization
clear all
% Working Directory
cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 2\Homework\Solution'

%% Question 1: Part 1

rng(2021)

% Consumers
I = 2000;

mu = 1;
sigma = 1;
alpha = lognrnd(mu,sigma^2,[I,1]); 

%exp_alpha = exp(mu+0.5*sigma^2);
%var_alpha = (exp(sigma^2)-1)*exp(2*mu+sigma^2);

%%% Firms
J = 4;
ownership = eye(J);

% consumer taste zeta_j net of price effect: 4*1
delta_np = zeros(J,1);
for j = 1:J
    delta_np(j) = j/5;
end

% marginal cost: 4*1
mc = zeros(J,1);
for j = 1:J
    mc(j) =  j/10;
end
 

%% Question 1: Part 2


opt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 10000,'MaxIter',5000);
grid = 30;
multi = unifrnd(0,2,grid,1); 


price_init = multi'.*mc;


price_noRT = zeros(J,grid);



for t = 1:grid
price_noRT(:,t) = fsolve(@(price)sys_foc(price, delta_np, alpha, I, ownership, mc), price_init(:,t), opt);
end

Var_price_grid = var(price_noRT,0,2);

sj = zeros(J,1);
profit = zeros(J,1);

[sj,profit] = profit_firms(delta_np,price_noRT(:,1),alpha,I,mc);

%% Question 2: Part 1 Exogenous Outside Option

ownership_RT = ones(J);

price_RTExo = zeros(J, grid);

for t = 1:grid
price_RTExo(:,t) = fsolve(@(price)sys_foc(price, delta_np, alpha, I,ownership_RT , mc), price_init(:,t),  opt);
end

Var_price_RT_grid = var(price_RTExo,0,2);

sj_RTExo = zeros(J,1);
profit_RTExo = zeros(J,1);

[sj_RTExo,profit_RTExo] = profit_firms(delta_np,price_RTExo(:,1),alpha,I,mc);

%% Question 2: Part 2 Endogeneous Outside Option



%% Export Results into Latex

Table_Simulation_noRT = zeros(J,3);

Table_Simulation_noRT(:,1) = price_noRT(:,1);
Table_Simulation_noRT(:,2) = sj(:,1);
Table_Simulation_noRT(:,3) = profit(:,1);

Firm = 1:4;
noRTResultVarNames = {'Price','Market Share','Profit'};
noRTResultRowNames = append("Product ",string(Firm))';

Table_Simulation_noRT   = array2table(Table_Simulation_noRT,'VariableNames',noRTResultVarNames,'RowNames',noRTResultRowNames);



Table_Simulation_RTExo = zeros(J,3);

Table_Simulation_RTExo(:,1) = price_RTExo(:,1);
Table_Simulation_RTExo(:,2) = sj_RTExo(:,1);
Table_Simulation_RTExo(:,3) = profit_RTExo(:,1);

RTExoResultRowNames = append("Product ",string(Firm));
RTExoResultVarNames = {'Price','Market Share','Profit'};
Table_Simulation_RTExo   = array2table(Table_Simulation_RTExo,'VariableNames',RTExoResultVarNames,'RowNames',RTExoResultRowNames);

table2latex(Table_Simulation_noRT,'Table_Simulation_noRT.tex');
table2latex(Table_Simulation_RTExo,'Table_Simulation_RTExo.tex');
