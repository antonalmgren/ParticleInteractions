clear all 
close all
%%
H = 50;
load('epsilon_7.mat')
epsilon7 = SA;
load('epsilon_5.mat')
epsilon5 = SA;
load('epsilon_4.mat')
epsilon4 =SA;
load('alpha_01.mat')
base = SA;
load('alpha_03.mat')
alpha03 = SA;
load('alpha_05.mat')
alpha05 = SA;
load('alpha_07.mat')
alpha07 = SA;
load('a_16.mat')
a16 = SA;
load('a_17.mat')
a17 = SA;
load('a_18.mat')
a18 = SA;
load('a_19.mat')
a19 = SA;
load('a_20.mat')
a20 = SA;
load('a_21.mat')
a21 = SA;
load('a_22.mat')
a22 = SA;
load('a_23.mat')
a23 = SA;

load('x4.mat');
load('x5.mat');

SAstring = [{'epsilon4'},{'epsilon5'},{'epsilon7'},{'base'},{'alpha03'},{'alpha05'},{'alpha07'},{'a18'},{'a19'},{'a20'},{'a22'},{'a23'}];
%% Export flux
nR = length(x4);
load('DELTA5.mat')
DELTA5 = DELTA;
load('DELTA4.mat')
DELTA4 = DELTA;


    
exportDT = epsilon4.wWhites.*epsilon4.M/H;
exportDT_x(1,:) = sum(exportDT,1)./DELTA;
exportFlux(1) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = epsilon5.wWhites.*epsilon5.M/H;
exportDT_x(2,:) = sum(exportDT,1)./DELTA;
exportFlux(2) = sum(exportDT,'all','omitnan')*1E-3; %[mgC/m2/d]

exportDT = epsilon7.wWhites.*epsilon7.M/H;
exportDT_x(3,:) = sum(exportDT,1)./DELTA;
exportFlux(3) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = base.wWhites.*base.M/H;
exportDT_x(4,:) = sum(exportDT,1)./DELTA;
exportFlux(4) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = alpha03.wWhites.*alpha03.M/H;
exportDT_x(5,:) = sum(exportDT,1)./DELTA;
exportFlux(5) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = alpha05.wWhites.*alpha05.M/H;
exportDT_x(6,:) = sum(exportDT,1)./DELTA;
exportFlux(6) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]


exportDT = alpha07.wWhites.*alpha07.M/H;
exportDT_x(7,:) = sum(exportDT,1)./DELTA;
exportFlux(7) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a16.wWhites.*a16.M/H;
exportDT_x(8,:) = sum(exportDT,1)./DELTA5;
exportFlux(8) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]


exportDT = a17.wWhites.*a17.M/H;
exportDT_x(9,:) = sum(exportDT,1)./DELTA5;
exportFlux(9) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a18.wWhites.*a18.M/H;
exportDT_x(10,:) = sum(exportDT,1)./DELTA5;
exportFlux(10) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a19.wWhites.*a19.M/H;
exportDT_x(11,:) = sum(exportDT,1)./DELTA;
exportFlux(11) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a20.wWhites.*a20.M/H;
exportDT_x(12,:) = sum(exportDT,1)./DELTA;
exportFlux(12) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a21.wWhites.*a21.M/H;
exportDT_x(13,:) = sum(exportDT,1)./DELTA;
exportFlux(13) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a22.wWhites.*a22.M/H;
exportDT_x(14,:) = sum(exportDT,1)./DELTA;
exportFlux(14) = sum(exportDT,'all','omitnan')*1E-3 ;%[mgC/m2/d]

exportDT = a23.wWhites.*a23.M/H;
exportDT_x(15,:) = sum(exportDT,1)./DELTA;
exportFlux(15) = sum(exportDT,'all','omitnan')*1E-3; %[mgC/m2/d]
%% Epsilon 
figure
semilogx(x4,exportDT_x(1:2,:),'LineWidth',1.5)
hold on
semilogx(x4,exportDT_x(4,:),'LineWidth',1.5)
semilogx(x4,exportDT_x(3,:),'LineWidth',1.5)
title('export flux')
text(0.01, 0.95, 'Export flux [mgC/m^2/d]' ,'Units','normalized')
text(0.01, 0.9,['\epsilon = 10^{-4}: ', num2str(round(exportFlux(1)))],'Units','normalized')
text(0.01, 0.85, ['\epsilon = 10^{-5}: ', num2str(round(exportFlux(2)))],'Units','normalized')
text(0.01, 0.8, ['\epsilon = 10^{-6}: ', num2str(round(exportFlux(4)))],'Units','normalized')
text(0.01, 0.75, ['\epsilon = 10^{-7}: ', num2str(round(exportFlux(3)))],'Units','normalized')
xlabel('radius [\mu m]')
ylabel('\mu g C m^{-2} d^{-1} \mu m^{-1}') 
legend('\epsilon=10^{-4}','\epsilon=10^{-5}','\epsilon=10^{-6}','\epsilon=10^{-7}','location','w')
set(gca,'FontSize',14)

%print('./plots/flux_epsilon.png', '-dpng', '-r400')
%% alpha


figure
semilogx(x4,exportDT_x(4:7,:),'LineWidth',1.5)
title('export flux')
text(0.01, 0.95, 'Export flux [mgC/m^2/d]' ,'Units','normalized')
text(0.01, 0.9,['\alpha = 0.1: ', num2str(round(exportFlux(4)))],'Units','normalized')
text(0.01, 0.85, ['\alpha = 0.3: ', num2str(round(exportFlux(5)))],'Units','normalized')
text(0.01, 0.8, ['\alpha = 0.5: ', num2str(round(exportFlux(6)))],'Units','normalized')
text(0.01, 0.75, ['\alpha = 0.7: ', num2str(round(exportFlux(7)))],'Units','normalized')
xlabel('radius [\mu m]')
ylabel('\mu g C m^{-2} d^{-1}\mu m^{-1}') 
legend('\alpha = 0.1','\alpha = 0.3','\alpha = 0.5','\alpha = 0.7','location','w')
set(gca,'FontSize',14)

%print('./plots/flux_alpha.png', '-dpng', '-r400')

%% a

%figure('Position', [300, 300, 700, 425])
semilogx(x5,exportDT_x(8:10,:),'LineWidth',1.5)
hold on
semilogx(x4,exportDT_x(11:15,:),'LineWidth',1.5)
title('export flux')
text(0.01, 0.95, 'Export flux [mgC/m^2/d]' ,'Units','normalized')
text(0.01, 0.9,['a = 1.6: ', num2str(round(exportFlux(8)))],'Units','normalized')
text(0.01, 0.9,['a = 1.7: ', num2str(round(exportFlux(8)))],'Units','normalized')
text(0.01, 0.85,['a = 1.8: ', num2str(round(exportFlux(9)))],'Units','normalized')
text(0.01, 0.8, ['a = 1.9: ', num2str(round(exportFlux(10)))],'Units','normalized')
text(0.01, 0.75, ['a = 2.0: ', num2str(round(exportFlux(11)))],'Units','normalized')
text(0.01, 0.7, ['a = 2.1: ', num2str(round(exportFlux(12)))],'Units','normalized')
text(0.01, 0.65, ['a = 2.2: ', num2str(round(exportFlux(13)))],'Units','normalized')
text(0.01, 0.6, ['a = 2.3: ', num2str(round(exportFlux(14)))],'Units','normalized')
xlabel('radius [\mu m]')
ylabel('\mu g C m^{-2} d^{-1}\mu m^{-1}') 
legend('a = 1.6','a = 1.7','a = 1.8','a = 1.9','a = 2.0','a = 2.1','a = 2.2','a = 2.3','location','sw')
set(gca,'FontSize',14)

%print('./plots/flux_a.png', '-dpng', '-r400')

%%
load('m.mat')
DELTA = DELTA*1E-4; %convert to cm

N = epsilon4.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3

NNN(1,:) = Nc./DELTA;
%%%%%
N = epsilon5.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(2,:) = Nc./DELTA;

N = base.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(3,:) = Nc./DELTA;

N = epsilon7.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(4,:) = Nc./DELTA;

N = base.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(5,:) = Nc./DELTA;

N = alpha03.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(6,:) = Nc./DELTA;

N = alpha05.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(7,:) = Nc./DELTA;

N = alpha07.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(8,:) = Nc./DELTA;

N = a16.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x5(2:end)-x5(1:end-1));
NNN(9,:) = Nc./DELTA;

N = a17.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x5(2:end)-x5(1:end-1));
NNN(10,:) = Nc./DELTA;

N = a18.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x5(2:end)-x5(1:end-1));
NNN(11,:) = Nc./DELTA;

N = a19.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(12,:) = Nc./DELTA;

N = a20.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(13,:) = Nc./DELTA;

N = a21.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(14,:) = Nc./DELTA;

N = a22.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(15,:) = Nc./DELTA;

N = a23.M./m;
Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA(1) = 1E-4*x4(1);
DELTA(2:nR) = 1E-4*(x4(2:end)-x4(1:end-1));
NNN(16,:) = Nc./DELTA;



slope  = @(x,b)   1E-4*x.^(-b); 


monterey = @(x) 0.027*x.^(-2.965);


%figure('Position', [300, 300, 700, 425])
loglog(2*x4*1E-4,NNN(1:4,:), 'LineWidth',1.5)
hold on
plot(2*x4*1E-4,slope(2*x4*1E-4,2))
plot(2*x4*1E-4,slope(2*x4*1E-4,4),':')
title('Particle size spectrum')
xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
legend('\epsilon=10^{-4}','\epsilon=10^{-5}','\epsilon=10^{-6}','\epsilon=10^{-7}','slope=-2','slope=-4','Location','SouthWest')
set(gca,'FontSize',16)

%print('./plots/spectrum_epsilon.png', '-dpng', '-r400')

%figure('Position', [300, 300, 700, 425])
loglog(2*x4*1E-4,NNN(5:8,:), 'LineWidth',1.5)
hold on
plot(2*x4*1E-4,slope(2*x4*1E-4,2))
plot(2*x4*1E-4,slope(2*x4*1E-4,4),':')
title('Particle size spectrum')
xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
legend('\alpha = 0.1','\alpha = 0.3','\alpha = 0.5','\alpha = 0.7','slope=-2','slope=-4','Location','SouthWest')
set(gca,'FontSize',16)

%print('./plots/spectrum_alpha.png', '-dpng', '-r400')

%figure('Position', [300, 300, 700, 425])
loglog(2*x5*1E-4,NNN(9:11,:), 'LineWidth',1.5)
hold on
loglog(2*x4*1E-4,NNN(12:16,:), 'LineWidth',1.5)
plot(2*x4*1E-4,slope(2*x4*1E-4,2))
plot(2*x4*1E-4,slope(2*x4*1E-4,4),':')

title('Particle size spectrum')
xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
legend('a = 1.6','a = 1.7','a = 1.8','a = 1.9','a = 2.0','a = 2.1','a = 2.2','a = 2.3','slope=-2','slope=-4','Location','SouthWest')
set(gca,'FontSize',16)

%print('./plots/spectrum_a.png', '-dpng', '-r400')
%%

figure('Position', [0, 0, 900, 750])
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'tight';
nexttile

semilogx(x4*2*1E-4,exportDT_x(1:2,:),'LineWidth',1.5)
hold on
semilogx(x4*2*1E-4,exportDT_x(4,:),'LineWidth',1.5)
semilogx(x4*2*1E-4,exportDT_x(3,:),'LineWidth',1.5)
title('Export flux')
% text(0.01, 0.95, 'Export flux [mgC/m^2/d]' ,'Units','normalized')
% text(0.01, 0.9,['\epsilon = 10^{-4}: ', num2str(round(exportFlux(1)))],'Units','normalized')
% text(0.01, 0.85, ['\epsilon = 10^{-5}: ', num2str(round(exportFlux(2)))],'Units','normalized')
% text(0.01, 0.8, ['\epsilon = 10^{-6}: ', num2str(round(exportFlux(4)))],'Units','normalized')
% text(0.01, 0.75, ['\epsilon = 10^{-7}: ', num2str(round(exportFlux(3)))],'Units','normalized')
%xlabel('radius [\mu m]')
ylabel('\mu g C m^{-2} d^{-1} \mu m^{-1}') 
xlim([2*min(x5)*1E-4 2*max(x5)*1E-4])
xticks([1E-3 1E-2 1E-1 1])
%legend('\epsilon=10^{-4}','\epsilon=10^{-5}','\epsilon=10^{-6}','\epsilon=10^{-7}','location','northEast')
set(gca,'FontSize',12)

nexttile

loglog(2*x4*1E-4,NNN(1:4,:), 'LineWidth',1.5)
hold on
plot(2*x4*1E-4,slope(2*x4*1E-4,2),'b--','LineWidth',1.5)
plot(2*x4*1E-4,slope(2*x4*1E-4,4)*1000,'k--','LineWidth',1.5)
title('Particle size spectrum')
%xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
xlim([2*min(x5)*1E-4 2*max(x5)*1E-4])
xticks([1E-3 1E-2 1E-1 1])
legend('\epsilon=10^{-4}','\epsilon=10^{-5}','\epsilon=10^{-6}','\epsilon=10^{-7}','slope=-2','slope=-4','Location','EastOutside')
% set(gca,'FontSize',12, 'YaxisLocation','right')
set(gca,'FontSize',12)

nexttile

semilogx(x4*2*1E-4,exportDT_x(4:7,:),'LineWidth',1.5)
%title('export flux')
% text(0.01, 0.95, 'Export flux [mgC/m^2/d]' ,'Units','normalized')
% text(0.01, 0.9,['\alpha = 0.1: ', num2str(round(exportFlux(4)))],'Units','normalized')
% text(0.01, 0.85, ['\alpha = 0.3: ', num2str(round(exportFlux(5)))],'Units','normalized')
% text(0.01, 0.8, ['\alpha = 0.5: ', num2str(round(exportFlux(6)))],'Units','normalized')
% text(0.01, 0.75, ['\alpha = 0.7: ', num2str(round(exportFlux(7)))],'Units','normalized')
%xlabel('radius [\mu m]')
ylabel('\mu g C m^{-2} d^{-1}\mu m^{-1}') 
xlim([2*min(x5)*1E-4 2*max(x5)*1E-4])
xticks([1E-3 1E-2 1E-1 1])
%legend('\alpha = 0.1','\alpha = 0.3','\alpha = 0.5','\alpha = 0.7','location','NorthEast')
set(gca,'FontSize',12)

nexttile

loglog(2*x4*1E-4,NNN(5:8,:), 'LineWidth',1.5)
hold on
plot(2*x4*1E-4,slope(2*x4*1E-4,2),'b--','LineWidth',1.5)
plot(2*x4*1E-4,slope(2*x4*1E-4,4)*100,'k--','LineWidth',1.5)
%title('Particle size spectrum')
%xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
xlim([2*min(x5)*1E-4 2*max(x5)*1E-4])
xticks([1E-3 1E-2 1E-1 1])
legend('\alpha = 0.1','\alpha = 0.3','\alpha = 0.5','\alpha = 0.7','slope=-2','slope=-4','Location','EastOutside')
% set(gca,'FontSize',12, 'YaxisLocation','right')
set(gca,'FontSize',12)

nexttile

semilogx(x5*2*1E-4,exportDT_x(8:10,:),'LineWidth',1.5)
hold on
semilogx(x4*2*1E-4,exportDT_x(11:15,:),'LineWidth',1.5)
%title('export flux')
% text(0.01, 0.95, 'Export flux [mgC/m^2/d]' ,'Units','normalized')
% text(0.01, 0.9,['a = 1.6: ', num2str(round(exportFlux(8)))],'Units','normalized')
% text(0.01, 0.9,['a = 1.7: ', num2str(round(exportFlux(8)))],'Units','normalized')
% text(0.01, 0.85,['a = 1.8: ', num2str(round(exportFlux(9)))],'Units','normalized')
% text(0.01, 0.8, ['a = 1.9: ', num2str(round(exportFlux(10)))],'Units','normalized')
% text(0.01, 0.75, ['a = 2.0: ', num2str(round(exportFlux(11)))],'Units','normalized')
% text(0.01, 0.7, ['a = 2.1: ', num2str(round(exportFlux(12)))],'Units','normalized')
% text(0.01, 0.65, ['a = 2.2: ', num2str(round(exportFlux(13)))],'Units','normalized')
% text(0.01, 0.6, ['a = 2.3: ', num2str(round(exportFlux(14)))],'Units','normalized')
xlabel('Particle Diameter [cm]')
ylabel('\mu g C m^{-2} d^{-1}\mu m^{-1}') 
xlim([2*min(x5)*1E-4 2*max(x5)*1E-4])
xticks([1E-3 1E-2 1E-1 1])
%legend('a = 1.6','a = 1.7','a = 1.8','a = 1.9','a = 2.0','a = 2.1','a = 2.2','a = 2.3','location','NorthEast')
set(gca,'FontSize',12)

nexttile

loglog(2*x5*1E-4,NNN(9:11,:), 'LineWidth',1.5)
hold on
loglog(2*x4*1E-4,NNN(12:16,:), 'LineWidth',1.5)
plot(2*x4*1E-4,slope(2*x4*1E-4,2),'b--','LineWidth',1.5)
plot(2*x4*1E-4,slope(2*x4*1E-4,4)*100,'k--','LineWidth',1.5)
%title('Particle size spectrum')
xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
xlim([2*min(x5)*1E-4 2*max(x5)*1E-4])
xticks([1E-3 1E-2 1E-1 1])
legend('a = 1.6','a = 1.7','a = 1.8','a = 1.9','a = 2.0','a = 2.1','a = 2.2','a = 2.3','slope=-2','slope=-4','Location','EastOutside')
% set(gca,'FontSize',12, 'YaxisLocation','right')
set(gca,'FontSize',12)