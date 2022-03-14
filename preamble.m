function [p,f] = preamble
% Loading parameters and inline functions to coagulation.m

%% Constants
p.rMax = 1E4; %[\mu m] max radius
p.rMin = 1 ; %[\mu m] min radius
p.a = 1.9; %Self similarity parameter
rho_sw = 1.027E-6; % density of seawater [\mug \mu m^-3] (from andy)
nu = 1E-6; % [m^2 s^-1] kinematic viscosity of seawater (from andy)
mu = nu*rho_sw*10^9;% [kg m^-1 s^-1  ] absolute viscosity (10^9 is a conversion factor for rho to kg/m^3) 
p.alpha =0.1; %stickiness
kb = 1.38065E-23; %Boltzmann constant [m^2 kg s^-2 K^-1]
p.epsilon = 1E-6; % [m^2 s^-3] %energy dissipation rate (McCave 1984) (converted from 1E-4 cm^2 s^-3)(1E-8)
p.remin = 0.1; % [d^-1] Remineralisation rate

%% environmental variables
T = 281; %temperature
p.H = 50; %[m] depth of mixed layer

%% grid and combination variables
p.nR = 20; %number of size bins
p.nD = 10; %number of density bins

deltaR = exp(log(p.rMax)/(p.nR-1)); % size step
deltaRho = 0.1*rho_sw/(p.nD-1);% 0.6*rho_sw/nD; %density step

p.q = deltaR^(p.a-3); % is this still valid when a is not part of delta?

p.L = p.nR*p.nD; %number of bins: b index [0, 1, ..., L-1]
K = (p.L+1)*p.L/2; % number of combos: k index [0, 1, ..., K-1]
k = [0:K-1]';
%b = [0:L-1]';
z = (2*p.L + 1 - sqrt((2*p.L + 1).^2 - 8*k))/2;%(2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2; % SHOULD IT BE 2L or L? is this eq 4.12?
p.bi = floor(z);
p.bj =   k - p.bi*p.L + p.bi.*(p.bi - 1)/2 + p.bi;%k - bi*L + bi.*(bi-1)/2 + L (4.13)
xi = floor(p.bi/p.nD); zi = p.bi - xi*p.nD; % Does this make sense? M is nD and N is nR
xj = floor(p.bj/p.nD); zj = p.bj - xj*p.nD;
p.x = [0:p.nR-1]; p.z = [0:p.nD-1];
[p.xMesh,p.zMesh] = meshgrid(p.x,p.z);

%% transformations

p.pip = 4*pi/3;

f.xz = @(b) [floor(b/(p.nD)), b - floor(b/(p.nD))*(p.nD)]; % bin number into x-z coordinates
f.logd = @(x) log(x)./log(deltaR); %used for finding daughter particles
f.zeta = @(xi,xj) (1 + deltaR.^(p.a*(xi - xj)));
f.p = @(x,z) x*(p.nD) + z; %gives bin number in a vector

f.r = @(x) p.rMin*deltaR.^x; % radius from x ordinate
f.y = @(x,z) deltaRho*z.*p.q.^x; % rho-rho_sw from x and z ordinate
f.mass = @(x,z) p.pip*(f.y(x,z) + rho_sw).*(f.r(x).^3);

w_func = @(x,z) (2.*1E9*f.y(x,z).*9.81.*(f.r(x)*1E-6).^2)./(9*mu);% [m/s]   *24*3600; %[m d^-1] sinking velocity NB! from de la rocha and passow 2007, unit double checked

%coordinates of the daughter particles for every combination
xioj = xi + f.logd(f.zeta(xj,xi))/p.a;
zioj = zi./f.zeta(xj,xi) + zj./f.zeta(xi,xj);

%% Target bins
% Dividing mass between four target bins
x300 = floor(xioj); %lowest x ordinate
z300 = floor(zioj); %lowest z ordinate
p.b300 = f.p(x300,z300);    %defining all four target bins (bin number)
p.b310 = f.p(x300+1,z300);
p.b301 = f.p(x300,z300+1);
p.b311 = f.p(x300+1,z300+1);
dx1 = xioj - x300;    %dividing mass between x and z ordinates
dx0 = 1 - dx1;
dz1 = zioj - z300;
dz0 = 1 - dz1;
p.f00 = dx0.*dz0; %determining the fraction going into each bin
p.f10 = dx1.*dz0;
p.f01 = dx0.*dz1;
p.f11 = dx1.*dz1;
% Keep inside state space (boundary condition)
p.b300(p.b300>(p.L-1)) = p.L-1;
p.b301(p.b301>(p.L-1)) = p.L-1;
p.b310(p.b310>(p.L-1)) = p.b310(p.b310>(p.L-1))-p.nD;
p.b311(p.b311>(p.L-1)) = p.b311(p.b311>(p.L-1))-p.nD;
p.b311(p.b311>(p.L-1)) = p.b311(p.b311>(p.L-1)) -1;



%% derived properties
%W = w_func(p.xMesh,p.zMesh)*24*3600;%[m d^-1]
    
w_tmp = w_func(p.xMesh,p.zMesh);
Re = 2E-6.*f.r(p.xMesh).*w_tmp./nu; %Reynolds number
fRe = 24./Re + 6./(1+Re.^0.5)+0.4; % whites approximation
w_it = sqrt((8E-6*f.r(p.xMesh).*9.81*1E9.*f.y(p.xMesh,p.zMesh))./(3E9*rho_sw.*fRe));
tick = 0;
while max(w_tmp./w_it,[],'all')>1 %iterate until converged (within 0.5-2 times the previous)
    w_tmp =w_it;
    Re = 2E-6.*f.r(p.xMesh).*w_it./nu;
    fRe = 24./Re + 6./(1+Re.^0.5)+0.4;
    w_it = sqrt((8E-6*f.r(p.xMesh).*9.81*1E9.*f.y(p.xMesh,p.zMesh))./(3E9*rho_sw.*fRe));
    
    
    tick=tick+1;
end
p.wWhites = w_it*24*3600; %[m/d]



% Vector version for parent particles, used for differential settling
for i = 1:K
    wVeci(i,:) = w_it(zi(i)+1,xi(i)+1);
    wVecj(i,:) = w_it(zj(i)+1,xj(i)+1);
end

%% coagulation kernels
% Brownian motion
beta_b = (2*kb*T)./(3*mu)*(f.r(xi)+f.r(xj)).^2./(f.r(xi).*f.r(xj));  % [m^3/s], double checked

% Shear
pp = f.r(xi)./f.r(xj);
beta_s = 9.8*(pp.^2./(1+2*pp.^2)).*(p.epsilon/nu)^0.5.*(1E-6*(f.r(xi)+f.r(xj))).^3;  %[m^3/s], double checked

% differential settling
beta_d = 0.5*pi*(1E-6*f.r(xi)).^2.*abs(wVeci-wVecj); % [m^3/s], double checked
 
p.beta = (beta_b + beta_s + beta_d)*3600*24; %[m^3 d^-1]

%% Fragmentation


pfrag = linspace(0.001,0.5,p.nR);
pfrag = repmat(pfrag,p.nD,1);
pfrag = pfrag./(rho_sw+f.y(p.xMesh,p.zMesh))*1E-6;
denfrag = linspace(1,0.5,p.nD);
for i = 1:p.nR
    pfrag(:,i) = pfrag(:,i).*denfrag';
end
p.pfrag = pfrag(:);
p.frag_div = 0.5;


end
