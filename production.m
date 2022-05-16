function [prod,prod_tot] = production(M,k,H,rMin,deltaR,deltaRho,q,rho_sw,P)

% k=1 means case with diatoms and generalists
% k=2 is same as 1, but with no density difference
% k=3  is production over all density classes in the smallest size bin
% k=4 is production provided from external source (model/measurements)

prod_tot = 1E5; % Total production (0.1 g/m2/d) 

if k == 1 || k == 2 
    prod = zeros(size(M));

    prod(1:2,1) = 2*prod_tot;
    prod(2:3,2) = 2*prod_tot; 
    prod(3:4,3) = 2*prod_tot;
    prod(4:5,4) = 2*prod_tot;
    prod(5:6,5) = 2*prod_tot;
    
    prod(end-2:end,8:10) = 3*prod_tot;
    prod(end-1:end,11:12) = 3*prod_tot;
    prod(end,13:15) = 3*prod_tot;

    if k == 2
        tmp = prod;
        prod = zeros(size(M));
        prod(1:2,1) = 2*prod_tot;
        prod(2:3,2) = 2*prod_tot; 
        prod(3:4,3) = 2*prod_tot;
        prod(4:5,4) = 2*prod_tot;
        prod(5:6,5) = 2*prod_tot;
        prod(9,8:15) = sum(tmp(:,8:15),1);
    end

    prod = prod/H/(sum(prod,"all")/prod_tot);
elseif k == 3
    prod = zeros(size(M));
    prod(:,1) = prod_tot; 
    prod = prod/H/(sum(prod,"all")/prod_tot);
elseif k == 4
    prod = zeros(size(M));

    if height(P) ==2
        P(3,:) = 1088.4;
    end

    idxS = floor(log(P(2,:)/rMin)/log(deltaR))+1; %from r(x) function
    idxRho = floor((P(3,:)*1E-9-rho_sw)./(deltaRho*q.^(idxS-1)))+1; % from y(x,z) function

    idxS(idxS>size(M,2)) = size(M,2);
    idxRho(idxRho>size(M,1)) = size(M,1);

    for i = 1:width(P)    
        prod(idxRho(i),idxS(i)) = prod(idxRho(i),idxS(i))+P(1,i);

    end


    end
end