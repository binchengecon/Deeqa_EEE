function surplus = consumer(delta, sigma, p, nu, alpha)

ns = size(nu,2);
alternatives = repmat(delta+p*alpha,1,ns)+sigma*p*nu;
pricescaling = -(alpha+sigma*nu);
logsum = log(1+sum(exp(alternatives),1));


surplus = 1/ns*sum(logsum./pricescaling);


end