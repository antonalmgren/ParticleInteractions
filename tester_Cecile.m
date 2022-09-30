load ('sim_HL_noDensity.mat');

% 
% for i = 1:length(sim.Mtrans)
%     Mt(:,:,i) = reshape(sim.Mtrans(i,:),size(sim.M));
% end
% 
% 
%     
n= sum(sim.N,1);
r = sim.r;

M = sim.M;
m = sim.m;
m = sum(m.*M,1)./sum(M,1);

rho = sim.y+1.027E-6;
rho = sum(rho.*M,1)./sum(M,1);

w = sum(sim.w.*M,1)./sum(M,1);

export = w.*sum(M,1);
export_frac = export/(sum(sim.prod,'all')*sim.H);
figure
semilogx(r,export)
title('export')
figure
semilogx(r,export_frac)
title('Export frac')
%%
figure
loglog(r,sum(M,1))
%%

figure
loglog(r,sum(M,1)./m)
hold on
loglog(r,n);





