
clear all 
close all
% first part to be used in parameters.m
%% constants
rMax = 1E4; %[\mu m] max radius
rMin = 1 ; %[\mu m] min radius
a = 2; %fractal dimension
rho_sw = 1.027E-6; % density of seawater [\mug \mu m^-3] (from andy)
nu = 1E-6; % [m^2 s^-1] kinematic viscosity of seawater (from andy)
mu = nu*rho_sw*10^9;% [kg m^-1 s^-1  ] absolute viscosity (10^9 is a conversion factor for rho to kg/m^3) 
alpha =0.1; %stickiness
kb = 1.38065E-23; %Boltzmann constant [m^2 kg s^-2 K^-1]
epsilon = 1E-7; % [m^2 s^-3] %energy dissipation rate (McCave 1984) (converted from 1E-4 cm^2 s^-3)(1E-8)
remin = 0.1; % [d^-1] Remineralisation rate

H = 50; %[m] depth of mixed layer

%% grid and combination variables
nR = 20; %number of size bins
nD = 10; %number of density bins

deltaR = exp(log(rMax)/(nR-1)); % size step
%deltaR = 1.1;
%deltaR = 10^(log10(rMax)/nR);
%deltaR = 2^(1/a);
deltaRho = 0.1*rho_sw/(nD-1);% 0.6*rho_sw/nD; %density step

q = deltaR^(a-3); % is this still valid when a is not part of delta?

L = nR*nD; %number of bins: b index [0, 1, ..., L-1]
K = (L+1)*L/2; % number of combos: k index [0, 1, ..., K-1]
k = [0:K-1]';
b = [0:L-1]';
z = (2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;%(2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2; % SHOULD IT BE 2L or L? is this eq 4.12?
bi = floor(z);
bj =   k - bi*L + bi.*(bi - 1)/2 + bi;%k - bi*L + bi.*(bi-1)/2 + L (4.13)
xi = floor(bi/nD); zi = bi - xi*nD; % Does this make sense? M is nD and N is nR
xj = floor(bj/nD); zj = bj - xj*nD;
x = [0:nR-1]; z = [0:nD-1];
[xMesh,zMesh] = meshgrid(x,z);

%% transformations

pip = 4*pi/3;

xz = @(b) [floor(b/(nD)), b - floor(b/(nD))*(nD)]; % bin number into x-z coordinates
logd = @(x) log(x)./log(deltaR); %used for finding daughter particles
zeta = @(xi,xj) (1 + deltaR.^(a*(xi - xj)));
p = @(x,z) x*(nD) + z; %gives bin number in a vector

r = @(x) rMin*deltaR.^x; % radius from x ordinate
y = @(x,z) deltaRho*z.*q.^x; % rho-rho_sw from x and z ordinate
mass = @(x,z) pip*(y(x,z) + rho_sw).*(r(x).^3);

w_func = @(x,z) (2.*1E9*y(x,z).*9.81.*(r(x)*1E-6).^2)./(9*mu);% [m/s]   *24*3600; %[m d^-1] sinking velocity NB! from de la rocha and passow 2007, unit double checked

%coordinates of the daughter particles for every combination
xioj = xi + logd(zeta(xj,xi))/a;
zioj = zi./zeta(xj,xi) + zj./zeta(xi,xj);

%% Target bins
% Dividing mass between four target bins
x300 = floor(xioj); %lowest x ordinate
z300 = floor(zioj); %lowest z ordinate
b300 = p(x300,z300);    %defining all four target bins (bin number)
b310 = p(x300+1,z300);
b301 = p(x300,z300+1);
b311 = p(x300+1,z300+1);
dx1 = xioj - x300;    %dividing mass between x and z ordinates
dx0 = 1 - dx1;
dz1 = zioj - z300;
dz0 = 1 - dz1;
f00 = dx0.*dz0; %determining the fraction going into each bin
f10 = dx1.*dz0;
f01 = dx0.*dz1;
f11 = dx1.*dz1;
% Keep inside state space
b300(b300>(L-1)) = L-1;
b301(b301>(L-1)) = L-1;
b310(b310>(L-1)) = b310(b310>(L-1))-nD;
b311(b311>(L-1)) = b311(b311>(L-1))-nD;
b311(b311>(L-1)) = b311(b311>(L-1)) -1;

%% environmental variables
T = 281; %temperature

%% derived properties
%w = @(r,rho) (2.*1E9*(rho-rho_sw).*9.81.*(r*1E-6).^2)./(9*mu);% [m/s]   *24*3600; %[m d^-1] sinking velocity NB! from de la rocha and passow 2007, unit double checked
W = w_func(xMesh,zMesh)*24*3600;%[m d^-1]
    


w_tmp = w_func(xMesh,zMesh);
Re = 2E-6.*r(xMesh).*w_tmp./nu;
fRe = 24./Re + 6./(1+Re.^0.5)+0.4;
w_it = sqrt((8E-6*r(xMesh).*9.81*1E9.*y(xMesh,zMesh))./(3E9*rho_sw.*fRe));
tick = 0;
while max(w_tmp./w_it,[],'all')>1 %iterate until converged (within 0.5-2 times the previous)
    w_tmp =w_it;
    Re = 2E-6.*r(xMesh).*w_it./nu;
    fRe = 24./Re + 6./(1+Re.^0.5)+0.4;
    w_it = sqrt((8E-6*r(xMesh).*9.81*1E9.*y(xMesh,zMesh))./(3E9*rho_sw.*fRe));
    
    
    tick=tick+1;
end
wWhites = w_it*24*3600; %[m/d]



% Vector version for parent particles, used for differential settling
for i = 1:K
    wVeci(i,:) = w_it(zi(i)+1,xi(i)+1);
    wVecj(i,:) = w_it(zj(i)+1,xj(i)+1);
end





%% coagulation kernels
% Brownian motion
beta_b = (2*kb*T)./(3*mu)*(r(xi)+r(xj)).^2./(r(xi).*r(xj));  % [m^3/s], double checked

% Shear
pp = r(xi)./r(xj);
beta_s = 9.8*(pp.^2./(1+2*pp.^2)).*(epsilon/nu)^0.5.*(1E-6*(r(xi)+r(xj))).^3;  %[m^3/s], double checked

% differential settling
beta_d = 0.5*pi*(1E-6*r(xi)).^2.*abs(wVeci-wVecj); % [m^3/s], double checked
 
beta = (beta_b + beta_s + beta_d)*3600*24; %[m^3 d^-1]


%% beta rectilinear turbulent shear

%beta = 1.3*(epsilon/nu)^0.5*(1E-6*(r(xi)+r(xj))).^3*24*3600;
% %%
%for i = 1:nR
%    for j = 1:nR
%        betaplot(i,j)= 1.3*(epsilon/nu)^0.5*(1E-6*(r(x(i))+r(x(j)))).^3*24*3600;
%    end
%end
%% Fragmentation


pfrag = linspace(0.001,0.5,nR);
pfrag = repmat(pfrag,nD,1);
pfrag = pfrag./(rho_sw+y(xMesh,zMesh))*1E-6;
denfrag = linspace(1,0.5,nD);
for i = 1:nR
    pfrag(:,i) = pfrag(:,i).*denfrag';
end
pfrag = pfrag(:);
frag_div = 0.5;

%% Interactions
N = zeros(nD,nR); %number of particles/m^3 
M = zeros(nD,nR); % [\mug C/ m^3 
m = mass(xMesh,zMesh) ;
prod_tot = 1E5; %0.1 g/m2/d
prod = zeros(size(M));

prod(:,1) = prod_tot/nD;

% prod(1:3,1) = prod_tot/10/H; 
% prod(4:10,6) = prod_tot/10/H;
% prod(3,11)=prod_tot/10/H;
% 
% if nR == 20 && nD == 10
%     load('./init/M_20_10.mat')
%     disp('loaded M_20_10.mat')
% elseif nR ==30 && nD ==15
%     load('./init/M_30_15.mat')
%     disp('loaded M_30_15.mat')
% elseif length(beta)==1
%     load('./init/M_20_10_beta4.mat')
%     disp('loaded M_20_10_beta4.mat')
%     
% else
    M = prod;
    disp('without spinup')
% end
%load('./init/M_20_10_beta.mat')

N = M./m;

%% Transient solution 
tic
options = odeset('NonNegative',1:length(M(:)));
[t,dM] = ode23(@interactionsDT, [0:10000], [M(:) ],options,m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div);
runtime = toc        
M = reshape(dM(end,:),nD,nR);

RMSE = rms(dM(end,:)-dM(end-1,:))

SA.M = M;
SA.a = a;
SA.alpha = alpha;
SA.epsilon = epsilon;
%SA.beta = betaplot;
SA.RMSE = RMSE;
SA.wWhites = wWhites;

save('./init/M_20_10.mat','M')


%%

figure
surface(x,z,M)
title('M transient')
colorbar
set(gca,'ColorScale','log')

N = M./m;

figure
surface(x,z,N)
title('N transient')
colorbar
set(gca,'ColorScale','log')
set(gca,'ZScale','log')

for i = 1:length(t)
    MM(:,:,i) = reshape(dM(i,:),nD,nR);
end

exportDT = wWhites.*M/H;
exportDT_x = sum(exportDT,1);
exportFlux = sum(exportDT,'all','omitnan')*1E-3 %[mgC/m2/d]

figure
plot(x,exportDT_x)
title('export flux DT')
text(0.01, 0.95, ['Export flux ~ ',num2str(round(exportFlux)),'mgC/m^2/d'] ,'Units','normalized')
xlabel('size')
ylabel('\mu g C m^{-2} d^{-1}') 



%% for annual retreat

Nc = sum(N,1)*1E-6; % sum and change to #/cm^3
DELTA = 1E-4*r(1:nR);
NNN = Nc./DELTA;



slope  = @(x,b)   1E-4*x.^(-b); 


figure
loglog(2*r(x)*1E-4,NNN, 'LineWidth',2)
hold on
plot(2*r(x)*1E-4,slope(2*r(x)*1E-4,2))
plot(2*r(x)*1E-4,slope(2*r(x)*1E-4,3))
plot(2*r(x)*1E-4,slope(2*r(x)*1E-4,4))
title('Particle size spectrum')
xlabel('Particle Diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
legend('Size spectrum', 'slope = -2','slope = -3','slope = -4','Location','SouthWest')
set(gca,'FontSize',16)


%%
N3 = N;
N3(N3==0) = NaN;
N3 = log(N3);
%N3 = log10(N3);
N3 = N3+abs(min(N3,[],'all'));
N3(N3==0) = NaN;

w_weight = sum(wWhites.*N,1)./sum(N,1);

figure
loglog(r(x),w_weight,'-','LineWidth',2) 
hold on
scatter(r(xMesh(:)),wWhites(:),N3(:))
xlabel('Radius (\mu m)')
ylabel('Sinking velocity, (m/d)')
title('Mean sinking velocity of aggregates')


%Weight sinking speed by particle numbers


%% Diagnostics

for i = 1:length(t)
    Mt = reshape(dM(i,:),nD,nR);
    [Mdt,Mremin,Mfrag] = interactionsDT(t(i),Mt(:),m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div);
    frag(i) = sum(Mfrag);
    COM(i) = sum(Mdt)-sum(Mremin)+sum(Mt(:).*wWhites(:)/H)-sum(prod(:)) ;
end

figure
plot(t,COM)
%% w comparison

B117 = 1.10; % [cm^-0.17 s^-1] (Kries & Evans 2000)

B062 = (1.10^(-1/0.17))^(0.38);

w_PL = @(r,eta,B) B * (2*r.*1E-4).^eta*3600*24*1E-2;

W_PL117 = w_PL(r(x),1.17,B117);
W_PL062 = w_PL(r(x),0.62,B062);
%%
% 
% figure
% loglog(r(x),W_PL117,'r-',r(x),W_PL062,'g-',r(x),W,'b-',r(x),wWhites,'m-')
% % hold on
% % loglog(1E-3*r(1,:),W_PL117,'r-')
% %xlim([r(1,1) r(1,end)])
% %ylim([0 1000])
% legend('\eta = 1.17','\eta = 0.62', 'Stokes','Whites')
% ylabel('Sinking velocity [m d^{-1}]')
% xlabel('radius [\mu m]')

figure
loglog(r(x),W,'b-',r(x),wWhites,'r-')
 hold on
 loglog(r(x),W_PL117,'g-')
 loglog(r(x),W_PL062,'m-')
%xlim([1E2 r(x(end))])
%ylim([0 1000])
ylabel('Sinking velocity [m d^{-1}]')
xlabel('radius [\mu m]')
set(gca,'FontSize', 14)

%% Velocity figures
% 
% 
% 
% figure
% surface(Re)
% colorbar
% %set(gca,'xscale','log')
% set(gca,'ColorScale','log')
% title('Re')
% 
% figure
% subplot(1,3,1)
% surface(W)
% colorbar;
% title('w')
% %set(gca,'xscale','log')
% set(gca,'ColorScale','log')
% colormap jet
% 
% subplot(1,3,2)
% surface(wWhites)
% colorbar
% title('w Whites approximation')
% %set(gca,'xscale','log')
% set(gca,'ColorScale','log')
% colormap jet
% 
% subplot(1,3,3)
% surface(W./wWhites)
% colorbar
% title('w:w_whites')
% %set(gca,'xscale','log')
%  %set(gca,'ColorScale','log')
% colormap jet





