function [sim] = runCoagulation()
%%

[p,f]= preamble;
xz = f.xz;

N = zeros(p.nD,p.nR); %number of particles/m^3 
M = zeros(p.nD,p.nR); % [\mug C/ m^3 
p.m = f.mass(p.xMesh,p.zMesh) ;
prod_tot = 1E5; %0.1 g/m2/d
p.prod = zeros(size(M));

%prod(:,1) = prod_tot/nD;
p.prod(1,1) = 10*prod_tot;
 p.prod(2:3,1) = 2*prod_tot; 
 p.prod(4:10) = prod_tot; 
p.prod = p.prod/p.H/(sum(p.prod,"all")/prod_tot);
     M = p.prod;
    disp('without spinup')

    N = M./p.m;

    tic
options = odeset('NonNegative',1:length(M(:)));
[t,dM] = ode23(@interactions, [0:1000], [M(:) ],options,p,xz);
runtime = toc        
sim.M = reshape(dM(end,:),p.nD,p.nR);

sim.RMSE = rms(dM(end,:)-dM(end-1,:))

sim.p = p;

end
