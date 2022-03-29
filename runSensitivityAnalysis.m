clear all

a = [1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3];
alpha = [0.1,0.3,0.5,0.7];
epsilon = [1E-4,1E-5,1E-6,1E-7];
nR = 25;
nD = 10;
rMax = 1E4;
tMax = 50000;

for i = 1:length(a)

    [simA.sim(i)] = coagulation(a(i),0.1,1E-6,nR,nD,rMax,tMax);
end

for i = 1:length(alpha)
    [simAl.sim(i)] = coagulation(1.8,alpha(i),1E-6,nR,nD,rMax,tMax);
end

for i = 1:length(epsilon)
    [simE.sim(i)] = coagulation(1.8,0.1,epsilon(i),nR,nD,rMax,tMax);
end


save("sensitivityanalysis\a.mat","simA")
save("sensitivityanalysis\alpha.mat","simAl")
save("sensitivityanalysis\epsilon.mat","simE")

plotSensitivityAnalysis(simA,simAl,simE)