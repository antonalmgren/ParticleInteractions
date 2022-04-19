function [prod,prod_tot] = production(M,k,H)

% k=1 means case with diatoms and generalists
% k=2 means only 1 micron sized production

prod_tot = 1E5; % Total production (0.1 g/m2/d) 

if k ==1
    prod = zeros(size(M));

    prod(1,1) = prod_tot;
    prod(2:4,1) = prod_tot; 
    prod(8:10,6) = 3*prod_tot; 
    prod(9:10,7) = 3*prod_tot;
    prod(10,8) = 3*prod_tot;

    prod = prod/H/(sum(prod,"all")/prod_tot);
elseif k==2
    prod = zeros(size(M));
    prod(4,1:5) = prod_tot;
    prod = prod/H/(sum(prod,"all")/prod_tot);
end