clear all

a = [1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3];
alpha = [0.1,0.3,0.5,0.7];
epsilon = [1E-4,1E-5,1E-6,1E-7];
nR = 25;
nD = 10;
rMax = 1E5;

for i = 1:length(alpha)

    [aSim.sim(i)] = coagulation(a(i),0.1,1E-6,nR,nD,rMax);
end

% for i = 1:length(alpha)
%     [alphaSim.sim(i)] = coagulation(1.8,alpha(i),1E-6,nR,nD,rMax);
% end
% 
% for i = 1:length(epsilon)
%     [epsilonSim.sim(i)] = coagulation(1.8,0.1,epsilon(i),nR,nD,rMax);
% end


save("sensitivityanalysis\a.mat","aSim")
% save("sensitivityanalysis\alpha.mat","alphaSim")
% save("sensitivityanalysis\epsilon.mat","epsilonSim")

