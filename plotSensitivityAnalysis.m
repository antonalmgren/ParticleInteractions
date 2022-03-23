function  plotSensitivityAnalysis(simA,simAl,simE)

for i = 1:length(simA.sim)
    N = simA.sim(i).N;
    Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
    NspecA(i,:) = Nc./(simA.sim(i).DELTA*1E-4);
    exportA = (simA.sim(i).w.*simA.sim(i).M/simA.sim(i).H)/sum(simA.sim(i).prod,"all");
    exportAx(i,:) = sum(exportA,1,'omitnan')*100;
    labelsA(i) = {['a = ',num2str(simA.sim(i).a)]};
end



for i = 1:length(simAl.sim)
    N = simAl.sim(i).N;
    Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
    NspecAl(i,:) = Nc./(simAl.sim(i).DELTA*1E-4);
    exportAl = (simAl.sim(i).w.*simAl.sim(i).M/simAl.sim(i).H)/sum(simAl.sim(i).prod,"all");
    exportAlx(i,:) = sum(exportAl,1,'omitnan')*100;
    labelsAl(i) = {['\alpha = ',num2str(simAl.sim(i).alpha)]}; 
end

for i = 1:length(simE.sim)
    N = simE.sim(i).N;
    Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
    NspecE(i,:) = Nc./(simE.sim(i).DELTA*1E-4);
    exportE = (simE.sim(i).w.*simE.sim(i).M/simE.sim(i).H)/sum(simE.sim(i).prod,"all");
    exportEx(i,:) = sum(exportE,1,'omitnan')*100;
    labelsE(i) = {['\epsilon = ',num2str(simE.sim(i).epsilon)]}; 
end

slope  = @(x,b)   1E-5*x.^(-b); 

d = 2E-4*simA.sim(1).r;

figure('Position', [0, 0, 900, 750],'Color','white')
tiledlayout(3,2,TileSpacing="compact",Padding="tight")

nexttile
loglog(d,NspecA','LineWidth',1.5)
hold on
aSlope3 = loglog(d,slope(d,3),'k--','LineWidth',1.5);
aSlope4 = loglog(d,slope(d,4),'k:','LineWidth',1.5);
ylabel('Number spectrum [# cm^{-4}]')
%xlabel('Diameter [cm]')
legend([aSlope3, aSlope4],'slope = -3', 'slope = -4')
ylim([1E-2 1E10])
xlim([1E-4,1E-1])
title('A')

nexttile
semilogx(d,exportAx','LineWidth',1.5)
ylabel('Export [%]')
%xlabel('Diameter [cm]')
legend(labelsA)
xlim([1E-4,1E-1])
title('B')

nexttile
loglog(d,NspecAl','LineWidth',1.5)
hold on
alSlope3 = loglog(d,slope(d,3),'k--','LineWidth',1.5);
alSlope4 = loglog(d,slope(d,4),'k:','LineWidth',1.5);
ylabel('Number spectrum [# cm^{-4}]')
%xlabel('Diameter [cm]')
legend([alSlope3, alSlope4],'slope = -3', 'slope = -4')
ylim([1E-2 1E10])
xlim([1E-4,1E-1])
title('C')

nexttile
semilogx(d,exportAlx','LineWidth',1.5)
ylabel('Export [%]')
%xlabel('Diameter [cm]')
legend(labelsAl)
xlim([1E-4,1E-1])
title('D')

nexttile
loglog(d,NspecE','LineWidth',1.5)
hold on
eSlope3 = loglog(d,slope(d,3),'k--','LineWidth',1.5);
eSlope4 = loglog(d,slope(d,4),'k:','LineWidth',1.5);
ylabel('Number spectrum [# cm^{-4}]')
xlabel('Diameter [cm]')
legend([eSlope3, eSlope4],'slope = -3', 'slope = -4')
ylim([1E-2 1E10])
xlim([1E-4,1E-1])
title('E')

nexttile
semilogx(d,exportEx','LineWidth',1.5)
ylabel('Export [%]')
xlabel('Diameter [cm]')
legend(labelsE)
xlim([1E-4,1E-1])
title('F')

print('./figures/SensitivityAnalysis.png', '-dpng', '-r400')

end