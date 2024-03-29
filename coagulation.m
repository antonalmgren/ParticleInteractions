function [sim] = coagulation(a,alpha,epsilon,nR,nD,rMax,tMax,prodCase,P,M)
arguments
    a (1,1) {mustBeInRange(a,0,3)} = 1.8; %self similarity parameter
    alpha (1,1) {mustBeInRange(alpha,0,1)} =0.1; %stickiness
    epsilon (1,1) {mustBeInRange(epsilon,1E-9,1E-2)} = 1E-6; % [m^2 s^-3] %energy dissipation rate 
    nR (1,1) = 25; %number of size bins
    nD (1,1)= 10; %number of density bins
    rMax (1,1) = 1E4; %[\mu m] max radius
    tMax (1,1) = 10000; % No. of timesteps
    prodCase (1,1) = 4; % production structure, set in production.m
    P (:,:) = [ 1064E3,1302E3 ,1912E3 ,1859E3 ,4904E3 ;1,1,1,6.41,7.92;1088.4,1088.4,1088.4,1233.3,1100]; %production [\mug C/ m^2 /d ]
    M (:,:) = zeros(nD,nR); %initial condition [\mug C/ m^3 ]
end

    simulation_time = [0:tMax];


%% ------------------------------------------------------------------------
% constants
% -------------------------------------------------------------------------
rMin = 1 ; %[\mu m] min radius
rho_sw = 1.027E-6; % density of seawater [\mug \mu m^-3] (from andy)
nu = 1E-6; % [m^2 s^-1] kinematic viscosity of seawater (from andy)
mu = nu*rho_sw*10^9;% [kg m^-1 s^-1  ] absolute viscosity (10^9 is a conversion factor for rho to kg/m^3) 
kb = 1.38065E-23; %Boltzmann constant [m^2 kg s^-2 K^-1]
remin = 0.1; % [d^-1] Remineralisation rate

% Environmental parameters
H = 50; %[m] depth of mixed layer
T = 281; %temperature

%% ------------------------------------------------------------------------
% Grid and combination variables
% -------------------------------------------------------------------------
deltaR = exp(log(rMax)/(nR-1)); % size step
deltaRho = 2*rho_sw/(nD-1);% 0.6*rho_sw/nD; %density step

q = deltaR^(a-3); % is this still valid when a is not part of delta?

L = nR*nD; %number of bins: b index [0, 1, ..., L-1]
K = (L+1)*L/2; % number of combos: k index [0, 1, ..., K-1]
k = [0:K-1]';
b = [0:L-1]';
z = (2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;%(2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2; % SHOULD IT BE 2L or L? is this eq 4.12?
bi = floor(z);
bj =   k - bi*L + bi.*(bi - 1)/2 + bi;%k - bi*L + bi.*(bi-1)/2 + L (4.13)
xi = floor(bi/nD); zi = bi - xi*nD; 
xj = floor(bj/nD); zj = bj - xj*nD;
x = [0:nR-1]; z = [0:nD-1];
[xMesh,zMesh] = meshgrid(x,z);

%% ------------------------------------------------------------------------
% Transformations
% -------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Target bins
% -------------------------------------------------------------------------

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
while any(b310>L-1)
    b310(b310>(L-1)) = b310(b310>(L-1))-nD;
end
while any(b311>L)
    b311(b311>(L-1)) = b311(b311>(L-1))-nD;
end
b311(b311>(L-1)) = b311(b311>(L-1)) -1;



%% ------------------------------------------------------------------------
% Sinking Velocity, White's approximation
% -------------------------------------------------------------------------

W = w_func(xMesh,zMesh)*24*3600;%[m d^-1]

w_tmp = w_func(xMesh,zMesh); % Velocity using Stoke's law
Re = 2E-6.*r(xMesh).*w_tmp./nu; % Reynolds number
fRe = 24./Re + 6./(1+Re.^0.5)+0.4; % Drag correction for high Re
w_it = sqrt((8E-6*r(xMesh).*9.81*1E9.*y(xMesh,zMesh))./(3E9*rho_sw.*fRe)); %White's
tick = 0;

% Iterate until converged sinking velocity
while max(w_tmp./w_it,[],'all')>1 
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

%% ------------------------------------------------------------------------
% Coagulation kernels
% -------------------------------------------------------------------------

% Brownian motion
beta_b = (2*kb*T)./(3*mu)*(r(xi)+r(xj)).^2./(r(xi).*r(xj));  % [m^3/s], double checked

% Turbulent shear
pp = r(xi)./r(xj);
beta_s = 9.8*(pp.^2./(1+2*pp.^2)).*(epsilon/nu)^0.5.*(1E-6*(r(xi)+r(xj))).^3;  %[m^3/s], double checked

% differential settling
beta_d = 0.5*pi*(1E-6*r(xi)).^2.*abs(wVeci-wVecj); % [m^3/s], double checked
 
beta = (beta_b + beta_s + beta_d)*3600*24; %[m^3 d^-1]

%% ------------------------------------------------------------------------
% Fragmentation (simple version)
% -------------------------------------------------------------------------
pfrag = linspace(0.001,0.5,nR);
pfrag = repmat(pfrag,nD,1);
pfrag = pfrag./(rho_sw+y(xMesh,zMesh))*1E-6;
denfrag = linspace(1,0.5,nD);
for i = 1:nR
    pfrag(:,i) = pfrag(:,i).*denfrag';
end
pfrag = pfrag(:);
frag_div = 0.5;

%% ------------------------------------------------------------------------
% State variables and production
% -------------------------------------------------------------------------
N = zeros(nD,nR); %number of particles/m^3 
%M = zeros(nD,nR); % [\mug C/ m^3 ]
m = mass(xMesh,zMesh) ;

[prod,prod_tot] = production(M,prodCase,H,rMin,deltaR,deltaRho,q,rho_sw,P);


if sum(M,"all")==0
    M = prod;
end


%% ------------------------------------------------------------------------
% Transient solution 
% -------------------------------------------------------------------------

tic
disp('starting simulation')
options = odeset('NonNegative',1:length(M(:)));
[t,dM] = ode23(@interactions, simulation_time, [M(:) ],options,m,xMesh, ...
    zMesh,bi,bj,nR ,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha, ...
    beta,wWhites,L,H,prod,remin,pfrag,frag_div);
runtime = toc;
disp(['Simulation finished in ',num2str(runtime),' seconds'])

%% ------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------

sim.M = reshape(dM(end,:),nD,nR);
sim.N = sim.M./m;
sim.RMSE = rms(dM(end,:)-dM(end-1,:));
sim.a = a;
sim.alpha = alpha;
sim.epsilon = epsilon;
sim.nR = nR;
sim.nD = nD;
sim.H = H;
sim.m = m;
sim.w = wWhites;
sim.prod = prod;
sim.prod_tot = prod_tot;
sim.r = r(x);
sim.y = y(xMesh,zMesh);
sim.DELTA(1:nR-1) = r(x(2:nR))-r(x(1:nR-1));
sim.DELTA(nR) = r(nR+1)-r(nR);
sim.Mtrans = dM;
sim.remin = remin;
sim.export = sim.M.*sim.w/sim.H;

%% ------------------------------------------------------------------------
% Diagnostics
% -------------------------------------------------------------------------

for i = 1:length(t)
    Mt = reshape(dM(i,:),nD,nR);
    [Mdt,Mremin,Mfrag] = interactions(t(i),Mt(:),m,xMesh,zMesh,bi,bj,nR ...
        ,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites, ...
        L,H,prod,remin,pfrag,frag_div);
    frag(i) = sum(Mfrag);
    COM(i) = sum(Mdt)-sum(Mremin)+sum(Mt(:).*wWhites(:)/H)-sum(prod(:)) ;
end

figure
plot(t,COM)
hold on 
plot(t,frag)
xlabel('time step')
ylabel('error')
legend('COM','fragmentation loss')
end


