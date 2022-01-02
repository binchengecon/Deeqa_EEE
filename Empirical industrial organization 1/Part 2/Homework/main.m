%%
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
    mc(j) = j/10;
end
 


range = 10;
k = 1+ 0.1*(1:range); %try initial prices of 1.1 to 2.0 times marginal cost.
price_init = k.*mc;



opt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 10000,'MaxIter',5000);

price1 = zeros(J, size(k,2));
for t = 1:size(k,2)
price1(:,t) = fsolve(@(price)sys_foc(price, delta_np, alpha, I, ownership, mc), price_init(:,t), opt);
end
price1 % same values for all different initial prices.

%%
up = exp(delta_np + alpha'.*price_init(:,1)); % 5 x 2000
down = 1+sum(up,1); % 1 x 2000, col sum
s_hat_ij = up./down; %  5 x 2000

second = alpha'.*s_hat_ij*s_hat_ij';
first = diag(sum(alpha'.*s_hat_ij,2)); % row sum

sj = (1/I)*sum(s_hat_ij,2);
omega = (1/I)*(first-second);  % J X J matrix

res = 1e10*(sj + omega.*ownership * (price_init(:,1)-mc));

sys_foc(price_init(:,1), delta_np, alpha, I, ownership, mc)

%% Q2

price2 = zeros(J, size(k,2));

for t = 1:size(k,2)
price2(:,t) = fsolve(@(price)sys_foc2(price, delta_np, alpha, I, mc), price_init(:,t), opt);
end
price2

%summary
all_price = [price1(:,2) price2(:,2)]; %prices from initial guess of 1.2*mc
all_price

profit1 = profit(delta_np,all_price(:,1),alpha,I,mc);

profit2 = profit(delta_np,all_price(:,2),alpha,I,mc);
%% Q3

cs1 = consumer(delta_np,alpha,all_price(:,1)); %this is the total consumer surplus
cs2 = consumer(delta_np,alpha,all_price(:,2));
profit1 = profit(delta_np,all_price(:,1),alpha,I,mc); % total profit
profit2 = profit(delta_np,all_price(:,2),alpha,I,mc);

totalwelfare1=cs1+profit1;
totalwelfare2=cs2+profit2;

totalwelfare2-totalwelfare1
percent(totalwelfare2,totalwelfare1)

