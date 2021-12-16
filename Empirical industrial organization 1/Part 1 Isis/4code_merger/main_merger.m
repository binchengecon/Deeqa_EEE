%set seed so we all generate the same data
rng(7)

%define number of products
J = 4;

%generate a identifying number for each firm: 1,...4
ident_firm = [1:J]';

%draw costs in a U(2,3)
cost = 2 + rand(J,1);

%draw demand shocks in N(0,1)
xi = randn(J,1);

%set price coefficient
alpha = -0.8;

%compute mean utilities of product net of the price
%recall that x = intercept and beta = 1
delta_np = 1 + xi;


%% market outcomes under competition

%define ownership matrix under competition = identity
ownership = eye(J);

%set fsolve options
opt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 10000,'MaxIter',5000);

%define initial values for solving equilibrium prices 
price_init = 1.2*cost;

p_comp = fsolve(@(price)sys_foc(price, delta_np, alpha, cost, ownership), price_init, opt);

%compute market shares given the equilibrium prices
edelta = exp(delta_np + alpha*p_comp);

sum_edelta = 1+sum(edelta); %adds the outside good

sj = edelta./sum_edelta;

%market share of the outside option, defined by 1-sum(sj) or by 1./sum_edelta
s0 = 1-sum(sj);

%average mark up: 2 definitions
mean((p_comp-cost)./p_comp)

mean((p_comp)./cost)


%% market outcomes under the merger

%replace firm's identifier by one if it is equal to 2
ident_firm(ident_firm == 2) = 1;

%update ownership matrix
ownership_merger = (repmat(ident_firm, 1, J) == repmat(ident_firm', J, 1));

%cost efficiencies: merging parties benefit from synergies
cost_merger = cost;
cost_merger(1:2) = cost(1:2) - 0.8;

%solve for new equilibrium prices: same function applied to new cost and
%new ownership matrix
p_merger = fsolve(@(price)sys_foc(price, delta_np, alpha, cost_merger, ...
    ownership_merger), price_init, opt);

%average mark up
mean((p_merger-cost_merger)./p_merger)

mean(p_merger./cost_merger)



%% market outcome after the merger under partial collusion

%dummy colluding firms: 1-2 and 3 collude, firm 4 competes
ident_coll_firm = [1;1;1;0];

kappa = 0.8;

%matrix with the degree of internalization of other's products = 1 if firm
%= 1-2 and product = 3 or firm =3 and products = 1 or 2
internaliz =  (repmat(ident_coll_firm, 1, J) == repmat(ident_coll_firm', J, 1)) - ownership_merger;

%define new "ownership" matrix under partial collusion
ownership_coll = ownership_merger +  kappa * internaliz;

%solve new price after merger with partial collusion
p_merger_coll = fsolve(@(price)sys_foc(price, delta_np, alpha, cost_merger,...
    ownership_coll), price_init, opt);


price_time = [p_comp, p_merger, p_merger_coll];

time = [1,2,3];

figure
plot(time, price_time, 'LineWidth', 2)
legend('Product 1', 'Product 2', 'Product 3', 'Product 4')
