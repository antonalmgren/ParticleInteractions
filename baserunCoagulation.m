function[sim]= baserunCoagulation

clear all
close all

a = 1.8; %self similarity parameter
    alpha  =0.1; %stickiness
    epsilon = 1E-6; % [m^2 s^-3] %energy dissipation rate 
    nR  = 30; %number of size bins
    nD = 20; %number of density bins
    rMax = 1E5; %[\mu m] max radius
    tMax = 100; % No. of timesteps
    prodCase = 1; % production structure, set in production.m
    P  = 0; %production [\mug C/ m^2 /d ]
    
    %M  = zeros(nD,nR); %initial condition [\mug C/ m^3 ]
    load('mHLinit.mat');
prodCase = 4;
%P = [ 479E3,972E3 ,1544E3 ,1500E3 ,2200E3;1,1,1,6.41,7.92;1088.4,1088.4,1088.4,1233.3,1100]; %low light (van de poll 2013)
P = [ 1064E3,1302E3 ,1912E3 ,1859E3 ,4904E3 ;1,1,1,6.41,7.92;1088.4,1088.4,1088.4,1233.3,1100]; %High light
sim = coagulation(a,alpha,epsilon,nR,nD,rMax,tMax,prodCase,P,M)

plot0D(sim)
end