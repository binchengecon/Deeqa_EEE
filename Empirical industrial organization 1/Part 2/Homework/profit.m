function [a, output] = profit(delta_np,price,alpha,I,mc)

[~, sj]= mega(delta_np,price,alpha,I);

a = (price-mc).*sj*2000;
output = sum(a);
end

