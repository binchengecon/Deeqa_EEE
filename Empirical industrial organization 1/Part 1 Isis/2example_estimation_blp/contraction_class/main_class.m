%This code provides an example of the BLP contraction mapping
%we use fake data composed of a main dataset "data.csv" and the instruments
%"instru"
%The objective is to code the BLP market share inversion and the GMM
%function to estimate a random coefficient model
%note: the true sigma is 0.4

%% first: load the data and define the variables
data = csvread('data.csv');

z = csvread('instru.csv');

%market share vector
sj = data(:,1);

%market share of the outside option
s0 = data(:,2);

%exogenous product characteristics
xo_temp = data(:,3:7);

%price
price = data(:,8);

%index for the market
num_market = data(:,10);

%number of observations
kt = size(sj,1);

%add intercept to the exogenous characteristics
xo = [ones(kt,1), xo_temp];

%calcuate (z'z)^-1
invZ = (z'*z)\eye(size(z,2));

%%Useful indexes
%number of markets
T = max(num_market);

%vector with the number of products in each market
prodsMarket = zeros(T,1);
for tt = 1:T
    prodsMarket(tt) = sum(num_market==tt);
end

%marketStarts & marketEnds: locate beginning and end of each market
marketStarts = zeros(T,1);
marketEnds = zeros(T,1);
marketStarts(1) = 1;
marketEnds(1) = prodsMarket(1);

for tt = 2:T
    marketStarts(tt) = marketEnds(tt-1) + 1;
    marketEnds(tt) = marketStarts(tt) + prodsMarket(tt) - 1;
end

%compute an initial value for delta = delta under the simple logit model
delta_init = log(sj./s0);


%% Inversion of the market share
%draw 300 normal random draws to approximate the integral
ns = 300;

%draw in the normal(0, 1)
nu = randn(1, ns);

%try the contraction at the true sigma (=0.4)
%be careful, the function works for an inversion market by market so we
%need to loop over market
delta_new = zeros(kt,1);
for tt = 1: T
    Jt = marketStarts(tt):marketEnds(tt);
    delta_new(Jt) = blp_contraction(0.4, sj(Jt), price(Jt), delta_init(Jt), nu);
end

%% Estimation of the RC model

opt = optimset('display', 'iter', 'tolX', 1e-4);

%try to evaluate the GMM objective function at the value 0.1
gmm_init = gmm_blp(0.1, sj, xo, price, z, invZ, nu, delta_init, ...
    marketStarts, marketEnds);

gmm_fun = zeros(11,1);
jj = 0;
for ii = -1:0.1:1
    jj = jj + 1;
    gmm_fun(jj) = gmm_blp(ii, sj, xo, price, z, invZ, nu, delta_init, ...
        marketStarts, marketEnds);
end

%plot the GMM objective function between 0 and 1
figure
plot(-1:0.1:1, gmm_fun, 'k', 'linewidth', 1.7)
hold on
plot([0.4;0.4], [min(gmm_fun); max(gmm_fun)], 'r', 'linewidth', 1.7)
xlabel('\sigma')
ylabel('GMM objective function')
axis tight

%minimize the objective function to get sigma_hat
sigma_hat = fminunc(@(theta)gmm_blp(theta, sj, xo, price, z, invZ, nu, delta_init, ...
    marketStarts, marketEnds), 0.1, opt);







