function cs = consumer(delta_np,alpha,price)



a = exp(delta_np + alpha'.*price); % 5 x 2000
b = sum(a,1); %sum of products for each person i, colSum  1 x 2000


c = log(1+b); % add outside option 
d = c./abs(alpha'); %CS for each i
cs = sum(d);
