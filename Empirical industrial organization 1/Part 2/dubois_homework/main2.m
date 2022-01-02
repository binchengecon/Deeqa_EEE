clear all
clc

cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 2\dubois_homework'
%% Q1
ns = 200000;
%draw in the normal(0.1)
mu = 1;
sigma = 1;
alpha = lognrnd(mu,sigma^2,[ns,1])'; 

opt = optimset('display', 'iter', 'tolX', 1e-4);
%sigma

consumers = 2000;
%pure delta (without price effect)
delta_pure = [(4+(1/5)),(4+(2/5)),(4+(3/5)),(4+(4/5)),(5)]';
%price_int = [5,5,5,5,5]';
%marginal cost 
mc = [(1+(1/8)),(1+(2/8)),(1+(3/8)),(1+(4/8)),(1+(5/8))]';
price_int = 1.1*mc;
%onwership matrix
own_matrix = eye(5);
% price
price_Q1 = fsolve(@(price)help_foc(delta_pure, alpha, price, ns, own_matrix, mc), price_int, opt);
cs_Q1 = consumer_surplus(delta_pure,alpha,price_Q1,ns); % multiply by consumers ?
[no_need, sj_Q1] = omega_share(delta_pure,alpha,price_Q1,ns);
totalprofit_Q1 = sum((price_Q1-mc).*sj_Q1)*consumers;
wealfare_Q1 = cs_Q1*consumers + totalprofit_Q1;
%% Q2

%onwership matrix
own_matrix_new = [1,1,0,0,0;1,1,0,0,0;0,0,1,1,1;0,0,1,1,1;0,0,1,1,1];
% price
price_Q2 = fsolve(@(price)help_RPM_foc(delta_pure, alpha, price, ns, mc), price_int, opt);
cs_Q2 = consumer_surplus(delta_pure,alpha,price_Q2,ns);  % multiply by consumers ?
[no_need, sj_Q2] = omega_share(delta_pure,alpha,price_Q2,ns);
totalprofit_Q2 = sum((price_Q2-mc).*sj_Q2)*consumers;
wealfare_Q2 = cs_Q2*consumers + totalprofit_Q2;

%% Q2 endogenous
%%prices, market shares are the same as above
%%we need to compute counterfactual market share
retail_dummy = [1,1,1,1,1];
firstTerm = retail_dummy*((price_Q2-mc).*sj_Q2*consumers);

%%
firstTerm = firstTerm.*retail_dummy;
secondTerm = zeros(5,1);
for s=1:5
    delta_pure_cf = delta_pure;
    delta_pure_cf(s) = [];
    price_cf = price_Q2;
    price_cf(s) = [];
    mc_cf = mc;
    mc_cf(s)=[];
    [~,sj_cf_temp] = omega_share(delta_pure_cf,alpha,price_cf,ns);
    if s == 1
        sj_cf = [0; sj_cf_temp];
    elseif s== 5
        sj_cf = [sj_cf_temp; 0];
    else
        sj_cf = [sj_cf_temp(1:(s-1)) ;0; sj_cf_temp(s:end)];
    end
    secondTerm(s) = retail_dummy(:,s)' * retail_dummy*((price_Q2-mc).*sj_cf*consumers);
end
secondTerm = secondTerm'.*retail_dummy;
allFees = firstTerm - secondTerm;
manu_profit = sum(allFees,1);
retail_profit = retail_dummy*((price_Q2-mc).*sj_Q2*consumers)-sum(allFees,2);
%Check
totalprofit_Q2 - sum(retail_profit) - sum(manu_profit)
%% Q3
change_in_wealfare = wealfare_Q2 - wealfare_Q1;
%% results
price_Q1
price_Q2
change_in_wealfare
