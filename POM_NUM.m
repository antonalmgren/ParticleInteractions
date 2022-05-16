clear all
close all
%%
load("sim_LL_noDensity.mat")

%% Calculations

% stoke's sinking
rMesh = repmat(sim.r,sim.nD,1);
rho_sw = 1.027E-6; % density of seawater [\mug \mu m^-3] (from andy)
nu = 1E-6; % [m^2 s^-1] kinematic viscosity of seawater (from andy)
mu = nu*rho_sw*10^9;% [kg m^-1 s^-1  ] absolute viscosity (10^9 is a conversion factor for rho to kg/m^3) 
a = sim.a;
deltaR = exp(log(sim.r(end))/(sim.nR-1)); 
q = deltaR^(a-3);
rho_plankton = 1.0884E-6; % [\mug \mu m^-3]
x = 0:sim.nR-1;
y = (rho_plankton-rho_sw).*q.^x;


Fz = @ (z,Fz0,zRemin,H) Fz0.*exp((H-z)./zRemin);
z = 50:4000;

wR = @(r,y) (2.*1E9*y.*9.81.*(r*1E-6).^2)./(9*mu)*24*3600;

%%
Export = sim.export;
ExportSum = sum(Export,1,'omitnan')/sum(sim.prod,"all")*100;
ExportTot = sum(ExportSum)


f1  = fit(log10(sim.r)',ExportSum','gauss1');
%%
idx=0;
for i = 1:sim.nR
    FzFit(:,i) = Fz(z,f1(log10(sim.r(i))),wR(sim.r(i),y(i))/sim.remin,sim.H);
    for j = 1:sim.nD
        idx = idx+1;
        Fz1(:,idx) = Fz(z,Export(j,i),sim.w(j,i)/sim.remin,sim.H);

    end
end

Fz1Rel=sum(Fz1,2,"omitnan")/sum(sim.prod,'all')*100;
FzFitSum = sum(FzFit,2);

%%
figure
plot(Fz1Rel,z)
hold on
plot(FzFitSum,z)
axis ij

%%
ft = fittype('a*x^b'); 
N1d = sum(sim.N,1);

f2 = fit(sim.r',N1d',ft);


