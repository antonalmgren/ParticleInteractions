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
%%


%% Calculations

% stoke's sinking
rMesh = repmat(simLL.r,simLL.nD,1);
rho_sw = 1.027E-6; % density of seawater [\mug \mu m^-3] (from andy)
nu = 1E-6; % [m^2 s^-1] kinematic viscosity of seawater (from andy)
mu = nu*rho_sw*10^9;% [kg m^-1 s^-1  ] absolute viscosity (10^9 is a conversion factor for rho to kg/m^3) 
zRemin = @(sim) sim.w./sim.remin;
zRemin1 = zRemin(simLL);
Fz = @ (z,Fz0,zRemin,H) Fz0.*exp((H-z)./zRemin);
z = 50:4000;

rho_plankton = 1.0884E-6; % [\mug \mu m^-3]
wR = @(r) (2.*1E9*(rho_plankton-rho_sw).*9.81.*(r*1E-6).^2)./(9*mu)*24*3600;
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


f1  = fit(log10(simLLnD.r)',ExportSumLLnD','gauss1');

%% flux attenuation (from cavan 2019)
idx=0;
for i = 1:simLL.nR
    FzFit(:,i) = Fz(z,f1(log10(simLLnD.r(i))),wR(simLLnD.r(i))/simLLnD.remin,simLLnD.H);
    for j = 1:simLL.nD
        idx = idx+1;
        Fz1(:,idx) = Fz(z,ExportWhitesLL(j,i),zRemin1(j,i),simLL.H);
        Fz2(:,idx) = Fz(z,ExportWhitesLLnD(j,i),zRemin1(j,i),simLLnD.H);
        Fz3(:,idx) = Fz(z,ExportWhitesHL(j,i),zRemin1(j,i),simHL.H);
        FzMartinLL(:,idx) = ExportWhitesLL(j,i)*(z/simLL.H).^-0.87;
        Fz4(:,idx) = Fz(z,ExportWhitesHLnD(j,i),zRemin1(j,i),simHLnD.H);
        FzMartinHL(:,idx) = ExportWhitesHL(j,i)*(z/simHL.H).^-0.87;
        
    end
end


%Use 2d export flux to get the curve for each bin.

%% Plots



tiledlayout(3,2,'Padding',"tight","TileSpacing","compact")
nexttile([1 2])
semilogx(simLL.r*2E-4,ExportSumLL)
hold on
semilogx(simLLnD.r*2E-4,ExportSumLLnD)
semilogx(simHL.r*2E-4,ExportSumHL)
semilogx(simHLnD.r*2E-4,ExportSumHLnD)
legend('Low Light','Low Light 1D','High Light','High light 1D')
title('Export')

nexttile
imagesc(simLL.prod)
title('LL')
axis xy

nexttile
imagesc(simLLnD.prod)
title('LLnD')
axis xy

nexttile
imagesc(simHL.prod)
title('HL')
axis xy

nexttile
imagesc(simHLnD.prod)
title('HLnD')
axis xy

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

Fz1Rel=sum(Fz1,2,"omitnan")/sum(simLL.prod,'all')*100;

Fz2Rel=sum(Fz2,2,"omitnan")/sum(simLLnD.prod,'all')*100;
Fz3Rel=sum(Fz3,2,"omitnan")/sum(simHL.prod,'all')*100;
Fz4Rel=sum(Fz4,2,"omitnan")/sum(simHLnD.prod,'all')*100;
FzMartinRel=sum(FzMartinLL,2,"omitnan")/sum(simLL.prod,'all')*100;
FzMartinRel4=sum(FzMartinHL,2,"omitnan")/sum(simHL.prod,'all')*100;



figure
plot(Fz1Rel,z)
hold on
plot(Fz2Rel,z)
plot(Fz3Rel,z)
plot(Fz4Rel,z)
plot(FzMartinRel,z)
plot(FzMartinRel4,z)
plot(sum(FzFit,2),z)
axis ij
legend('LL','LLnD','HL','HLnD','martin LL','martin HL')
