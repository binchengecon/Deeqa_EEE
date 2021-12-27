%% Q1

%%% Consumers

I = 2000;
mu = 1;
sigma = 1;

rng(2020)
alpha = -lognrnd(mu,sigma^2,[I,1]); %1 X 2000, we have utility equ of +alpha NOT -alpha

%exp_alpha = exp(mu+0.5*sigma^2);
%var_alpha = (exp(sigma^2)-1)*exp(2*mu+sigma^2);

%%% Firms
J = 5;
ownership = eye(J);

% taste without price effect
delta_np = zeros(J,1);
for j = 1:J
    delta_np(j) = 4 + j/5;
end

%marginal cost
mc = zeros(J,1);
for j = 1:J
    mc(j) = 1 + j/8;
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

%% Q2

price2 = zeros(J, size(k,2));

for t = 1:size(k,2)
price2(:,t) = fsolve(@(price)sys_foc2(price, delta_np, alpha, I, mc), price_init(:,t), opt);
end
price2

%summary
all_price = [price1(:,2) price2(:,2)]; %prices from initial guess of 1.2*mc
all_price

%% Q3

cs1 = consumer(delta_np,alpha,all_price(:,1)); %this is the total consumer surplus
cs2 = consumer(delta_np,alpha,all_price(:,2));
profit1 = profit(delta_np,all_price(:,1),alpha,I,mc); % total profit
profit2 = profit(delta_np,all_price(:,2),alpha,I,mc);

totalwelfare1=cs1+profit1;
totalwelfare2=cs2+profit2;

totalwelfare2-totalwelfare1
percent(totalwelfare2,totalwelfare1)

