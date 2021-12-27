%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Homework\HW1
%
% Auto-generated by MATLAB on 2020-10-20 21:02:24
clearvars
cd 'C:\Users\33678\Desktop\Deeqa_EEE\Empirical industrial organization 1\Part 1 Isis\Homework\HW1'
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 25);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["year", "Var2", "num_manuf", "Var4", "num_brand", "Var6", "price", "cylinder1", "weight", "horsepower", "fuel_cost", "nb_households", "sj", "i1b_poids_vide", "i1b_puis_fisc", "i1b_cylindree", "i1b_euro_km", "i2b_poids_vide", "i2b_puis_fisc", "i2b_cylindree", "i2b_euro_km", "i3b_poids_vide", "i3b_puis_fisc", "i3b_cylindree", "i3b_euro_km"];
opts.SelectedVariableNames = ["year", "num_manuf", "num_brand", "price", "cylinder1", "weight", "horsepower", "fuel_cost", "nb_households", "sj", "i1b_poids_vide", "i1b_puis_fisc", "i1b_cylindree", "i1b_euro_km", "i2b_poids_vide", "i2b_puis_fisc", "i2b_cylindree", "i2b_euro_km", "i3b_poids_vide", "i3b_puis_fisc", "i3b_cylindree", "i3b_euro_km"];
opts.VariableTypes = ["double", "string", "double", "string", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var4", "Var6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var4", "Var6"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("base_project_deeqa.csv", opts);

%% Convert to output type
year = tbl.year;
num_manuf = tbl.num_manuf;
num_brand = tbl.num_brand;
price = tbl.price;
cylinder1 = tbl.cylinder1;
weight = tbl.weight;
horsepower = tbl.horsepower;
fuel_cost = tbl.fuel_cost;
nb_households = tbl.nb_households;
sj = tbl.sj;
i1b_poids_vide = tbl.i1b_poids_vide;
i1b_puis_fisc = tbl.i1b_puis_fisc;
i1b_cylindree = tbl.i1b_cylindree;
i1b_euro_km = tbl.i1b_euro_km;
i2b_poids_vide = tbl.i2b_poids_vide;
i2b_puis_fisc = tbl.i2b_puis_fisc;
i2b_cylindree = tbl.i2b_cylindree;
i2b_euro_km = tbl.i2b_euro_km;
i3b_poids_vide = tbl.i3b_poids_vide;
i3b_puis_fisc = tbl.i3b_puis_fisc;
i3b_cylindree = tbl.i3b_cylindree;
i3b_euro_km = tbl.i3b_euro_km;


%% Clear temporary variables
clear opts

%%%%% DATA
kt = size(sj,1);
xo_temp = [cylinder1 weight horsepower fuel_cost];
xo = [ones(kt,1), xo_temp]; % include intercept
z = [xo i1b_poids_vide i1b_puis_fisc i1b_cylindree i1b_euro_km i2b_poids_vide i2b_puis_fisc i2b_cylindree i2b_euro_km i3b_poids_vide i3b_puis_fisc i3b_cylindree i3b_euro_km];
invZ = inv(z'*z);


%% Q1.1



%%%% compute the share of outside option
T = max(year)-min(year)+1; % 6 years
prodsYear = zeros(T,1);

for t = 1:T
prodsYear(t) = sum(year==(2002+t));
end

yearStarts = zeros(T,1);
yearEnds = zeros(T,1);
yearStarts(1) = 1;
yearEnds(1) = prodsYear(1);

for t=2:T
    yearStarts(t) = yearEnds(t-1) + 1;
    yearEnds(t) = yearStarts(t) + prodsYear(t) - 1;
end

s0 = zeros(kt,1);
qq = zeros(T,1);

for tt = 1:T
    Jt = yearStarts(tt):yearEnds(tt);
    qq(tt) = sum(sj(Jt));
    s0(Jt) = 1-qq(tt);
end
clear qq


%%%%%% Estimation of the RC model

% sd normal 100 draws
rng(2020)
ns = 100;
nu = randn(1,ns);

% inital guess of delta
delta_init = log(sj./s0);

% optimize using gmm
opt = optimset('display', 'iter', 'tolX', 1e-4);
sigma_hat = fminunc(@(sigma)gmm_blp(sigma, sj, xo, price, z, invZ, nu, delta_init, ...
    yearStarts, yearEnds), 0.7, opt); % 0.7 initial guess of sigma
sigma_hat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JUST CHECKING : graphical check of sigma
gmm_fun = zeros(11,1);
jj = 0;
for ii = 0:0.1:1
    jj = jj +1;
    gmm_fun(jj) = gmm_blp(ii, sj, xo, price, z, invZ, nu, delta_init, yearStarts, yearEnds);
end
figure
plot(0:0.1:1, gmm_fun) % gmm objective function is the lowest around sigma 0.67
% fminunc seems working fine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% Inversion of the market share using sigma_hat
delta_new = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_new(Jt), pb_cv(tt)] = blp_contraction(sigma_hat, sj(Jt), price(Jt), delta_init(Jt), nu);   
end

[beta1, ~] = ivregression(delta_new, xo, price, z);





%% Q1.2

%%%% Dummy variable 

% dummy variable, year
tdum = zeros(kt,T);
for tt = 1:T
    Jt = yearStarts(tt):yearEnds(tt);
    tdum(Jt,tt)=1;
end
tdum (:, 1) = [];


% dummy brand
bdum = zeros(kt, max(num_brand));
for t=1:kt
    bdum(t,num_brand(t))=1;
end
bdum (:, 40) = [];
    
% xo = [ones(kt,1), xo_temp]
xo2 = [xo, tdum, bdum]; %add dummy to xo
z2 = [xo2 i1b_poids_vide i1b_puis_fisc i1b_cylindree i1b_euro_km i2b_poids_vide i2b_puis_fisc i2b_cylindree i2b_euro_km i3b_poids_vide i3b_puis_fisc i3b_cylindree i3b_euro_km];
invZ2 = inv(z2'*z2);
rcond(xo2'*xo2);
rcond(z2'*z2);

%%%%%% Estimation of the RC model

% optimize using gmm
opt = optimset('display', 'iter', 'tolX', 1e-4);
sigma_hat2 = fminunc(@(sigma)gmm_blp(sigma, sj, xo2, price, z2, invZ2, nu, delta_init, ...
    yearStarts, yearEnds), 2, opt); % 2 initial guess of sigma
sigma_hat2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JUST CHECKING : graphical check of sigma
gmm_fun2 = zeros(size(0:0.1:2.5,2),1);
jj = 0;
for ii = 0:0.1:2.5
    jj = jj +1;
    gmm_fun2(jj) = gmm_blp(ii, sj, xo2, price, z2, invZ2, nu, delta_init, yearStarts, yearEnds);
end
figure
plot(0:0.1:2.5, gmm_fun2) % gmm objective function is the lowest around sigma 2.25
% fminunc seems working fine?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% Inversion of the market share using sigma_hat
delta_new2 = zeros(kt,1);
for tt = 1: T %for each time 2003-2008
    Jt = yearStarts(tt):yearEnds(tt);  
    [delta_new2(Jt), pb_cv2(tt)] = blp_contraction(sigma_hat2, sj(Jt), price(Jt), delta_init(Jt), nu);   
end

[bta, ~] = ivregression(delta_new2, xo2, price, z2);
 % intercept, cyliner, weight, horsepower, fuel cost, price
beta2 =[bta(1:5)' bta(end)]';

result = [beta1, beta2]




%% 2 counterfactual

%data for 2008 only
Index8 = yearStarts(6):yearEnds(6);
tbl8 = tbl(Index8,:);
num_manuf8 = num_manuf(Index8);
num_brand8 = num_brand(Index8);
kt8 = size(Index8,2);
s08 = s0(Index8);
sj8 = sj(Index8);
price8 = price(Index8);
delta8 = delta_new2(Index8);
alpha = beta2(6); % alpha from Q1.2

%delta8 without price effect
delta_np = delta8 - alpha*price8;

%onwership matrix
own = (repmat(num_manuf8, 1, kt8) == repmat(num_manuf8', kt8, 1));
[omega, sj88]= mega(delta_np,nu,price8,sigma_hat2,alpha,ns);
omega_own = omega.*own;


%marginal cost
mc = price8 + omega_own \ sj8;

%average mark up
markup = mean((price8-mc)./price8); 

%total profit
totalprofit = sum((price8-mc).*sj8*nb_households(5028));

%average consumer surplus + constant
cs = consumer(delta_np,alpha,price8,sigma_hat2,nu,ns);

%% Q1.1
%volkswagen = 21, BMW = 1  num_manuf8

num_manuf8_merger = num_manuf8;
num_manuf8_merger(num_manuf8_merger==21) = 1;

own_merger = (repmat(num_manuf8_merger, 1, kt8) == repmat(num_manuf8_merger', kt8, 1));


optt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 10000,'MaxIter',5000);

% slow slow
p_merger = fsolve(@(price)sys_foc(price, delta_np, nu, sigma_hat2, alpha, ns, own_merger, mc), price8, optt);
% market share after merger
[~, sj_merger] = mega(delta_np,nu,p_merger,sigma_hat2,alpha,ns);

%%% average price effect, increases!
mean(p_merger)- mean(price8)
percent(mean(p_merger),mean(price8))

%%% average mark up
markup_merger = mean((p_merger-mc)./p_merger);
percent(markup_merger,markup)
markup_merger-markup

%%% profit of volkswagen + BMW
%index for Vol + BMW
indexVB = [1:kt8]';
for t = 1:kt8
    if num_manuf8_merger(t) == 1
        indexVB(t)=indexVB(t);
    else
        indexVB(t)=0;
    end
end
indexVB(indexVB==0) =[];


profitVB_merger = sum((p_merger(indexVB)-mc(indexVB)).*sj_merger(indexVB)*nb_households(5028));
profitVB = sum((price8(indexVB)-mc(indexVB)).*sj8(indexVB)*nb_households(5028));
percent(profitVB_merger,profitVB)

%%% total market profit
totalprofit_merger = sum((p_merger-mc).*sj_merger*nb_households(5028));
percent(totalprofit_merger,totalprofit)

%%% average CS + constant
cs_merger = consumer(delta_np,alpha,p_merger,sigma_hat2,nu,ns);
percent(cs_merger,cs)
cs_merger-cs



%% Q1.2

% guess lower bound is 60% of current mc
grid = [0.982:0.001:0.983];
gridsize = size(grid,2);
mc_grid = repmat(mc, 1, gridsize);

opttt = optimset('Display', 'iter',...
    'TolFun',1e-4,'TolCon',1e-3,'TolX',1e-4,...
    'MaxFunEvals', 15000,'MaxIter',5000);

for t = 1 : gridsize
mc_grid(indexVB,t)=grid(t)*mc(indexVB);
end

p_merger_grid = zeros(kt8, gridsize);
for t = 1:gridsize
    p_merger_grid(:,t) = fsolve(@(price)sys_foc(price, delta_np, nu, sigma_hat2, alpha, ns, own_merger, mc_grid(:,t)), p_merger, opttt);
end

cs_merger_grid = zeros(gridsize,1);
for t= 1:gridsize
    cs_merger_grid(t) = consumer(delta_np,alpha,p_merger_grid(:,t),sigma_hat2,nu,ns);
end

cs_merger_grid
%which one is closer to cs? at 0.983 !! reduction by 1.7% point
abs(cs-cs_merger_grid(1)) > abs(cs-cs_merger_grid(2))




%% Q2.1
% RENAULT = 15 buy 30% of PSA = 14
buy = 0.3;

%ownership matrix
num_manuf8_acq = num_manuf8;
num_manuf8_acq(num_manuf8_merger==14) = 15;
int0 =  (repmat(num_manuf8_acq, 1, kt8) == repmat(num_manuf8_acq', kt8, 1)) - own;

int1 = int0;

for i = 1:kt8
    for j= 1:kt8
        if i <j
            int1(i,j)=0;
        end
    end
end

own_acq = own + buy*int1;

% price after REN buy 30% of PSA
p_acq = fsolve(@(price)sys_foc(price, delta_np, nu, sigma_hat2, alpha, ns, own_acq, mc), price8, optt);

% market share after acq
[~, sj_acq] = mega(delta_np,nu,p_acq,sigma_hat2,alpha,ns);

%%% average price effect, increases!
mean(p_acq)- mean(price8)
percent(mean(p_acq),mean(price8))

%%% average markup
markup_acq = mean((p_acq-mc)./p_acq);
percent(markup_acq,markup)
markup_acq - markup

%%% profit of REN, PSA
%index for RENAULT
indexREN = [1:kt8]';
for t = 1:kt8
    if num_manuf8(t) == 15
        indexREN(t)=indexREN(t);
    else
        indexREN(t)=0;
    end
end
indexREN(indexREN==0) =[];

%index for PSA
indexPSA = [1:kt8]';
for t = 1:kt8
    if num_manuf8(t) == 14
        indexPSA(t)=indexPSA(t);
    else
        indexPSA(t)=0;
    end
end
indexPSA(indexPSA==0) =[];

profitPSA_acq = (1-buy)*sum((p_acq(indexPSA)-mc(indexPSA)).*sj_acq(indexPSA)*nb_households(5028));
profitPSA = sum((price8(indexPSA)-mc(indexPSA)).*sj8(indexPSA)*nb_households(5028));
profitREN_acq = sum((p_acq(indexREN)-mc(indexREN)).*sj_acq(indexREN)*nb_households(5028))+ buy/(1-buy)*profitPSA_acq;
profitREN = sum((price8(indexREN)-mc(indexREN)).*sj8(indexREN)*nb_households(5028));

%$$ price effect on PSA, REN
percent(mean(p_acq(indexPSA)),mean(price8(indexPSA)))
mean(p_acq(indexPSA))-mean(price8(indexPSA))

percent(mean(p_acq(indexREN)),mean(price8(indexREN)))
mean(p_acq(indexREN))-mean(price8(indexREN))

%$$ markup effect on PSA, REN
makeup_psa = mean((price8(indexPSA)-mc(indexPSA))./price8(indexPSA));
makeup_psa_acq = mean((p_acq(indexPSA)-mc(indexPSA))./p_acq(indexPSA));
percent(makeup_psa_acq,makeup_psa)
makeup_psa_acq-makeup_psa

makeup_ren = mean((price8(indexREN)-mc(indexREN))./price8(indexREN));
makeup_ren_acq = mean((p_acq(indexREN)-mc(indexREN))./p_acq(indexREN));
percent(makeup_ren_acq,makeup_ren)
makeup_ren_acq-makeup_ren
%$$ profit effect on PSA, REN
percent(profitPSA_acq,profitPSA)
profitPSA_acq - profitPSA

percent(profitREN_acq,profitREN)
profitREN_acq - profitREN

%%% total industry profit
totalprofit_acq = sum((p_acq-mc).*sj_acq*nb_households(5028));
percent(totalprofit_acq,totalprofit)
totalprofit_acq-totalprofit
%%% average CS + constant
cs_acq = consumer(delta_np,alpha,p_acq,sigma_hat2,nu,ns);
cs_acq-cs
percent(cs_acq,cs)



%% Q2.2

%if PSA also owns 30% of REN
own2=own*1;
own2(indexREN,indexREN)=0.7;
own2(indexPSA,indexPSA)=0.7; 

own_acq2 = own2 + buy*int0;

% price after REN PSA both buy 30%
p_acq2 = fsolve(@(price)sys_foc(price, delta_np, nu, sigma_hat2, alpha, ns, own_acq2, mc), price8, optt);

% market share after REN PSA both buy 30%
[~, sj_acq2] = mega(delta_np,nu,p_acq2,sigma_hat2,alpha,ns);

% both REN, PSA buy 30%, willingness
profitPSA_acq2 = (1-buy)*sum((p_acq2(indexPSA)-mc(indexPSA)).*sj_acq2(indexPSA)*nb_households(5028))+buy*sum((p_acq2(indexREN)-mc(indexREN)).*sj_acq2(indexREN)*nb_households(5028));
profitPSA_acq2-profitPSA_acq

profitREN_acq2 = (1-buy)*sum((p_acq2(indexREN)-mc(indexREN)).*sj_acq2(indexREN)*nb_households(5028))+buy*sum((p_acq2(indexPSA)-mc(indexPSA)).*sj_acq2(indexPSA)*nb_households(5028));
profitREN_acq2-profitREN_acq

%%% only PSA buy 30%
int2 = int0;

for i = 1:kt8
    for j= 1:kt8
        if i > j
            int2(i,j)=0;
        end
    end
end

own_acq3 = own + buy*int2;


% price after only PSA buy 30%
p_acq3 = fsolve(@(price)sys_foc(price, delta_np, nu, sigma_hat2, alpha, ns, own_acq3, mc), price8, optt);

% market share after only PSA buy 30%
[~, sj_acq3] = mega(delta_np,nu,p_acq3,sigma_hat2,alpha,ns);

% only PSA buy 30%, willingness
profitPSA_acq3 = sum((p_acq3(indexPSA)-mc(indexPSA)).*sj_acq3(indexPSA)*nb_households(5028))+buy*sum((p_acq3(indexREN)-mc(indexREN)).*sj_acq3(indexREN)*nb_households(5028));
profitPSA_acq3-profitPSA
