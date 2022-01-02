function cs = consumer_surplus(delta_pure,alpha,price,ns)
%delta_pure,alpha,price,ns
surplus_ij = exp(delta_pure - alpha.*price);
surplus_i = sum(surplus_ij,1); %sum of products for each i
surplus_i_and_outside_option = log(1+surplus_i);
alphai = - alpha;
adjusted = surplus_i_and_outside_option./abs(alphai); %adjusted consumer surplus for each i
cs = sum(adjusted)/ns;