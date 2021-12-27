%% Initialization

clear all

cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 2\Homework'


%% Q1

%%% Consumers

I = 2000;
mu = 1;
sigma = 1;

rng(2021)
alpha = -lognrnd(mu,sigma^2,[I,1]); %1 X 2000, we have utility equ of +alpha NOT -alpha

%exp_alpha = exp(mu+0.5*sigma^2);
%var_alpha = (exp(sigma^2)-1)*exp(2*mu+sigma^2);

%%% Firms
J = 4;
ownership = eye(J);

% taste without price effect
delta_np = zeros(J,1);
for j = 1:J
    delta_np(j) = j/5;
end

%marginal cost
mc = zeros(J,1);
for j = 1:J
    mc(j) =  j/10;
end
 




opt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 10000,'MaxIter',5000);

price1 = zeros(J,1);

price1 = fsolve(@(price)sys_foc(price, delta_np, alpha, I, ownership, mc), mc, opt);

price1

[sj,profit] = profit_two(delta_np,price1,alpha,I,mc);

%% Q2

price2 = zeros(J, 1);


price2 = fsolve(@(price)sys_foc2(price, delta_np, alpha, I, mc), mc,  opt);


%summary
all_price = [price1 price2]; %prices from initial guess of 1.2*mc
all_price


[profit1,Sum_profit1] = profit(delta_np,all_price(:,1),alpha,I,mc); % total profit
[profit2,Sum_profit2] = profit(delta_np,all_price(:,2),alpha,I,mc);
