%This function provides the BLP contraction mapping to invert the market
%shares under a random coefficient model with random coefficient on price
%only
%This function is valid ONLY for inversion market by market
%It outputs delta(sigma) = s^-1(sigma), whether the alogrithm converged
%and the number of iterations

function [delta_new, pb_cv,  iter] = blp_contraction(sigma, sj, price, delta_init, nu)

J = size(delta_init,1);

%initialize the vectors of delta, delta_old amd delta_new must be different
%at the begining to start the while-loop
delta_old = zeros(J,1);

delta_new = delta_init;

%we add a maximum number of iterations (=1000), if it is reached the algorithm
%stops
%we need to set up a tolerance level (here 1e-7) to represents when
%delta_old = delta_new
iter = 0;

while max(abs(delta_new - delta_old))> 1e-7 && iter < 1000
    
    iter = iter + 1;%increment the number of iterations
    
    delta_old = delta_new;%update the value of delta_old
    
    delta_new = delta_old + log(sj) - log(shares(delta_old, sigma, price, nu));
    %compute delta_new using BLP contraction mapping
    
end

pb_cv = 0;
if iter == 1000
    pb_cv = 1;
end


% disp(['no iter =', num2str(iter)])