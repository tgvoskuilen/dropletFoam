clear all
close all
clc

% Load and compare reaction mechanisms from Chemkin files

Ns = 41;   %select 25, 29, or 41

thermofile = strcat('chemkin',filesep,'therm',num2str(Ns,'%2d'),'.dat');

Tth = 200:5:2500;
species = ReadCHEMKINThermo(thermofile,Tth);
Ns =  length(fieldnames(species));


T=280;
dHTest = 4*interp1(Tth,species.HNO3.h,T)  ...
    - 4*interp1(Tth,species.NO2.h,T) ...
    - 2*interp1(Tth,species.H2O.h,T) - interp1(Tth,species.N2.h,T);


RFNA.rho = 1550;


MMH.rho = 870;

rhoT = MMH.rho + RFNA.rho;

P = 101325; %Pa

T = 280; %K
R = 8.314;
A = 1e3;
B = 0;
Ta = 500; %K

dG = interp1(Tth,species.CH3NHNH2.g,T) + 3*interp1(Tth,species.HNO3.g,T) ...
    - interp1(Tth,species.CH3ONO2.g,T) - 2*interp1(Tth,species.HONO.g,T) ...
    - 2*interp1(Tth,species.H2O.g,T) - interp1(Tth,species.N2.g,T);
Kp = exp(-dG/(R*T));
Kc = Kp*(P/(R*T))^6;

%       MMH   RFNA  CH3ONO2 H2O HONO  N2
Wi   = [46.07 57    77      18  47    28]; %kg/kmol
rhoi = [MMH.rho*0.35  RFNA.rho*0.65    0        0    0     0 ]; %kg/m3
xi = rhoi./Wi; %kmol/m3

kf = A*exp(-Ta/T);
kr = kf/Kc;
rf = kf*xi(1)*xi(2)^3;
rr = kr*xi(3)*xi(4)^2*xi(5)^2*xi(6);



