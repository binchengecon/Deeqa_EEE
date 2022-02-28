function pf = profit(n,z,fc,comp)

pf_active = 2*1.05^(2*(z-1)) * ( comp /(n + 1)^2 + (1-comp) /(1 + 1)^2 ) - fc;
pf= max(0,pf_active);

end

