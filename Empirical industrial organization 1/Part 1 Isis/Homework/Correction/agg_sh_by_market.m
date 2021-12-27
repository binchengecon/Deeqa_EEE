function [sj, sij] = agg_sh_by_market(delta, mu)

ns = size(mu, 2); %no of draws
kt = size(delta,1); %no of products

numer = exp(repmat(delta,1,ns) + mu); 


%only for market per market inversion
sum1 = 1 + sum(numer);

denom1 = repmat(sum1,kt,1); 


sij = numer./(denom1);   % simulated shares for each draw

sj= mean(sij,2);  