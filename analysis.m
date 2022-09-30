clear all
close all
%%
    a  = 2; %self similarity parameter
    alpha =0.1; %stickiness
    epsilon  = 1E-6; % [m^2 s^-3] %energy dissipation rate 
    nR = 30; %number of size bins
    nD = 20; %number of density bins
    rMax  = 1E4; %[\mu m] max radius
    tMax  = 10000; % No. of timesteps
    H = 50;


%% Loading
% sim1 = coagulation(a,alpha,epsilon,nR,nD,rMax,tMax,1);
% sim2 = coagulation(a,alpha,epsilon,nR,nD,rMax,tMax,2);
% sim3 = coagulation(a,alpha,epsilon,nR,nD,rMax,tMax,3);
% sim4 = load("sim_LL.mat");
load("sim_LL.mat")
simLL = sim;
load("sim_LL_noDensity.mat")
simLLnD = sim;
load("sim_HL.mat")
simHL = sim;
load("sim_HL_noDensity.mat")
simHLnD = sim;
load("AntonColormap.mat")

%%


%% Calculations

zRemin = @(sim) sim.w./sim.remin;
zRemin1 = zRemin(simLL);
Fz = @ (z,Fz0,zRemin,H) Fz0.*exp((H-z)./zRemin);
z = 50:4000;

Nc = sum(simLL.N,1); % #/m^3
NspecLL = Nc./(simLL.DELTA); % #/m^3/\mu m

Nc = sum(simLLnD.N,1); % #/m^3
NspecLLnD = Nc./(simLLnD.DELTA); % #/m^3/\mu m

Nc = sum(simHL.N,1); % #/m^3
NspecHL = Nc./(simHL.DELTA); % #/m^3/\mu m

Nc = sum(simHLnD.N,1); % #/m^3
NspecHLnD = Nc./(simHLnD.DELTA); % #/m^3/\mu m

slope  = @(x,b)  1E11* x.^(-b);

%% Export flux

ExportWhitesLL = simLL.export;
ExportSumLL = sum(ExportWhitesLL,1,'omitnan')/sum(simLL.prod,"all")*100;
ExportTotLL = sum(ExportSumLL)

ExportWhitesLLnD = simLLnD.export;
ExportSumLLnD = sum(ExportWhitesLLnD,1,'omitnan')/sum(simLLnD.prod,"all")*100;
ExportTotLLnD = sum(ExportSumLLnD)

ExportWhitesHL = simHL.export;
ExportSumHL = sum(ExportWhitesHL,1,'omitnan')/sum(simHL.prod,"all")*100;
ExportTotHL = sum(ExportSumHL)

ExportWhitesHLnD = simHLnD.M.*simHLnD.w/simHLnD.H;
ExportSumHLnD = sum(ExportWhitesHLnD,1,'omitnan')/sum(simHLnD.prod,"all")*100;
ExportTotHLnD = sum(ExportSumHLnD)



%% flux attenuation (from cavan 2019)
idx=0;
for i = 1:simLL.nR
    for j = 1:simLL.nD
        idx = idx+1;
        Fz1(:,idx) = Fz(z,ExportWhitesLL(j,i),zRemin1(j,i),simLL.H);
        Fz2(:,idx) = Fz(z,ExportWhitesLLnD(j,i),zRemin1(j,i),simLLnD.H);
        Fz3(:,idx) = Fz(z,ExportWhitesHL(j,i),zRemin1(j,i),simHL.H);
        FzMartinLL(:,idx) = ExportWhitesLL(j,i)*(z/simLL.H).^-0.87;
        Fz4(:,idx) = Fz(z,ExportWhitesHLnD(j,i),zRemin1(j,i),simHLnD.H);
        FzMartinHL(:,idx) = ExportWhitesHL(j,i)*(z/simHL.H).^-0.9695;
        
    end
end
Fz3Sum = sum(Fz3,2)/sum(simHL.prod,"all");
z1000 = 100:1000;
Martin100 = Fz3Sum(51)*(z1000/100).^-.9695 ;

%Use 2d export flux to get the curve for each bin.

%% Plots

% figure 2

figure
% tiledlayout(2,1)
% nexttile
surface(simHL.r,0:simHL.nD-1,simHL.M)
shading flat
%title('Mass')
c = colorbar;
set(gca,'Colorscale','log','xscale','log','TickDir','both')
%colormap("winter")
colormap(AntonColormap)
c.Label.String = '\mu g m^{-3}';
c.Label.FontSize = 12;
xlabel('radius [\mu m]',FontSize=12)
ylabel('z',FontSize=12)
ylim([0 simHL.nD-1])

saveas(gcf,'figures/stateSpace','png')
% nexttile
% surface(simHL.r,1:simHL.nD,simHL.M)
% shading flat
% title('Mass')
% c = colorbar;
% set(gca,'Colorscale','log','xscale','log')
% c.Label.String = '\mu g m^{-3}';

%% figure 3
%ticks= logspace(0, 5, 6);
%ticks = [10 1000 100000];
ticks = [1 100 10000];
figure("Position",[50,50,600,800])
tiledlayout(4,2,'Padding','tight',"TileSpacing","tight")


h(1) = nexttile;
surface( simLL.r,0:simLL.nD-1,simLL.prod)
set(gca,'xscale','log','TickDir','out')
title('LL')
%xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
%xticklabels([])
ylabel('z',FontSize=12)
ylim([0 simLL.nD-1])
shading flat
c = colorbar;
%colormap("winter")
c.Label.String = '\mu g m^{-2} d^{-1}';
caxis([0 5E6])
c.Label.FontSize = 12;
colormap(AntonColormap)
colorbar("off")

h(2) = nexttile;
surface( simLLnD.r,0:simLLnD.nD-1,simLLnD.prod)
set(gca,'xscale','log','TickDir','out')
title('LL 1D')
%xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
%xticklabels([])
%ylabel('z',FontSize=12)
yticklabels([])
ylim([0 simLL.nD-1])
shading flat
colormap(AntonColormap)
c = colorbar;
c.Label.String = '\mu g m^{-2} d^{-1}';
c.Label.FontSize = 12;
caxis([0 5E6])
%colormap("winter")


h(3) = nexttile;
surface( simHL.r,0:simHL.nD-1,simHL.prod)
set(gca,'xscale','log','TickDir','out')
title('HL')
xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
% xticklabels([])
ylabel('z',FontSize=12)
ylim([0 simHL.nD-1])
shading flat
c = colorbar;
c.Label.String = '\mu g m^{-2} d^{-1}';
c.Label.FontSize = 12;
caxis([0 5E6])
%colormap("winter")
colormap(AntonColormap)
colorbar("off")


h(4) = nexttile;
surface( simHLnD.r,0:simHLnD.nD-1,simHLnD.prod)
set(gca,'xscale','log','TickDir','out')
title('HL 1D')
xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
% xticklabels([])
%ylabel('z',FontSize=12)
yticklabels([])
ylim([0 simHLnD.nD-1])
shading flat
c = colorbar;
c.Label.String = '\mu g m^{-2} d^{-1}';
c.Label.FontSize = 12;
caxis([0 5E6])
%colormap("winter")
colormap(AntonColormap)

nexttile([1 2])
loglog(simLL.r,NspecLL,'LineWidth',1.5)
hold on
loglog(simLLnD.r,NspecLLnD,'LineWidth',1.5)
loglog(simHL.r,NspecHL,'LineWidth',1.5)
loglog(simHLnD.r,NspecHLnD,'LineWidth',1.5)
loglog(simHL.r,slope(simHL.r,4),'--k')
xlabel('radius [\mu m]',FontSize=12)
ylabel('# m^{-3} \mu m^{-1}',FontSize=12)
legend('LL','LL 1D', 'HL', 'HL 1D')

nexttile([1,2])
semilogx(simLL.r,ExportSumLL,'LineWidth',1.5)
hold on
semilogx(simLLnD.r,ExportSumLLnD,'LineWidth',1.5)
semilogx(simHL.r,ExportSumHL,'LineWidth',1.5)
semilogx(simHLnD.r,ExportSumHLnD,'LineWidth',1.5)
xlabel('radius [\mu m]',FontSize=12)
ylabel('flux [% of PP]',FontSize=12)
legend('LL','LL 1D','HL','HL 1D',Location='best')
%title('Export')

saveas(gcf,'figures/results','png')

% figure
% semilogx(sum(Fz1,2,"omitnan"),z)
% hold on
% semilogx(sum(Fz2,2,"omitnan"),z)
% semilogx(sum(Fz3,2,"omitnan"),z)
% semilogx(sum(Fz4,2,"omitnan"),z)
% semilogx(sum(FzMartin1,2,"omitnan"),z)
% legend('1','2','3','4','martin')
% axis ij
% ylim([0 1000])
%%
figure("Position",[50,50,600,800])
tiledlayout(4,2,'Padding','tight',"TileSpacing","tight")


h(1) = nexttile;
surface( simLL.r,0:simLL.nD-1,simLL.prod)
set(gca,'xscale','log','TickDir','out')
%title('LL')
%xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
xticklabels({})
ylabel('z',FontSize=12)
ylim([0 simLL.nD-1])
shading flat


h(2) = nexttile;
surface( simLLnD.r,0:simLLnD.nD-1,simLLnD.prod)
set(gca,'xscale','log','TickDir','out')
%title('LL 1D')
%xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
xticklabels([])
%ylabel('z',FontSize=12)
yticklabels([])
ylim([0 simLL.nD-1])
shading flat



h(3) = nexttile;
surface( simHL.r,0:simHL.nD-1,simHL.prod)
set(gca,'xscale','log','TickDir','out')
%title('HL')
xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
% xticklabels([])
ylabel('z',FontSize=12)
ylim([0 simHL.nD-1])
shading flat



h(4) = nexttile;
surface( simHLnD.r,0:simHLnD.nD-1,simHLnD.prod)
set(gca,'xscale','log','TickDir','out')
%title('HL 1D')
xlabel('radius [\mu m]',FontSize=12)
xticks(ticks)
% xticklabels([])
%ylabel('z',FontSize=12)
yticklabels([])
ylim([0 simHLnD.nD-1])
shading flat


set(h,'Colormap',AntonColormap,'Clim',[0 5E6])
cbh = colorbar(h(1));
cbh.Location = 'layout';
cbh.Layout.Tile = 'North';
cbh.Label.String = '\mu g C m^{-2} d^{-1}';
cbh.Label.FontSize = 12;

nexttile([1 2])
loglog(simLL.r,NspecLL,'LineWidth',1.5)
hold on
loglog(simLLnD.r,NspecLLnD,'LineWidth',1.5)
loglog(simHL.r,NspecHL,'LineWidth',1.5)
loglog(simHLnD.r,NspecHLnD,'LineWidth',1.5)
loglog(simHL.r,slope(simHL.r,4),'--k')
xlabel('radius [\mu m]',FontSize=12)
ylabel('# m^{-3} \mu m^{-1}',FontSize=12)
legend('LL','LL 1D', 'HL', 'HL 1D')

nexttile([1,2])
semilogx(simLL.r,ExportSumLL,'LineWidth',1.5)
hold on
semilogx(simLLnD.r,ExportSumLLnD,'LineWidth',1.5)
semilogx(simHL.r,ExportSumHL,'LineWidth',1.5)
semilogx(simHLnD.r,ExportSumHLnD,'LineWidth',1.5)
xlabel('radius [\mu m]',FontSize=12)
ylabel('flux [% of PP]',FontSize=12)
legend('LL','LL 1D','HL','HL 1D',Location='best')
%title('Export')

%% figure 4

Fz1Rel=sum(Fz1,2,"omitnan")/sum(simLL.prod,'all')*100;

Fz2Rel=sum(Fz2,2,"omitnan")/sum(simLLnD.prod,'all')*100;
Fz3Rel=sum(Fz3,2,"omitnan")/sum(simHL.prod,'all')*100;
Fz4Rel=sum(Fz4,2,"omitnan")/sum(simHLnD.prod,'all')*100;
%FzMartinRel=sum(FzMartinLL,2,"omitnan")/sum(simLL.prod,'all')*100;


tick = 1;
for i = 1:simHL.nD:simHL.nD*simHL.nR
    Fz3RelStep(:,tick)=sum(Fz3(:,i:i+simHL.nD-1),2,"omitnan")/sum(simHL.prod,'all')*100;
    tick = tick + 1;

end
FzMartinRelHL=sum(FzMartinHL,2,"omitnan")/sum(simHL.prod,'all')*100;

figure
semilogx(cumsum(Fz3RelStep(:,1:end-1),2),z)
hold on
FzSum = semilogx(Fz3Rel,z,'k','LineWidth',2);
FzMartin = semilogx(FzMartinRelHL,z,'--k','LineWidth',1.5);
ylim([0 1000])
xlabel('flux [% of PP]')
ylabel('depth [m]')
axis ij
legends = [FzSum,FzMartin];
legend(legends,'Total carbon flux','Martin curve fit (b=-0.9695','Location','best')


%%




figure('Position',[100 100 500 800])
tiledlayout(3,1,"TileSpacing","tight",Padding="tight")
nexttile(1)
surface(simHL.r,0:simHL.nD-1,simHL.export/sum(simHL.prod,'all')*100)
xlabel('radius [\mu m]',FontSize=12)
ylabel('z',FontSize=12)
c = colorbar;
c.Label.String = '% of PP';
%colormap("winter")
colormap(AntonColormap)
set(gca,'xscale','log')
shading flat
ylim([0 simHL.nD-1])

nexttile([2 1])
plot(cumsum(Fz3RelStep(:,1:end-1),2),z)
hold on
FzSum = plot(Fz3Rel,z,'k','LineWidth',2);
FzMartin = plot(FzMartinRelHL,z,'--k','LineWidth',1.5);
%martin2 = semilogx(Martin100*100,100:1000,'-k',LineWidth=3);
ylim([0 1000])
xlim([1 100])
xlabel('flux [% of PP]',FontSize=12)
ylabel('depth [m]',FontSize=12)
axis ij
legends = [FzSum,FzMartin];
legend(legends,'Total flux','Martin curve fit','Location','best',FontSize=12) % (b=-0.9695')

saveas(gcf,'figures/export','png')

%% nexttile([2 3])
figure
plot(Fz1Rel,z)
hold on
plot(Fz2Rel,z)
plot(Fz3Rel,z)
plot(Fz4Rel,z)
plot(FzMartinRelHL,z)
%plot(sum(FzFit,2),z)
axis ij
legend('LL','LLnD','HL','HLnD','martin b = -0.9695','Location','best')
ylim([0 1000])

