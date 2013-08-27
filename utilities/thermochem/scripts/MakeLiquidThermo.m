clear all
close all
clc

%TODO: Use L(T) expression
% L(T) = Lb*((Tc-T)/(Tc-Tb))^0.38 to get Cp change

% Load thermodynamic database
db = ReadDb();

% MMH
dH = 37.21; %kJ/mol, Schmidt p 296
Tb = 365; %K
dCp = -45; %J/mol/K (based on hydrazine)

MakeLiquid('CH3NHNH2L',db.CH3NHNH2, dH, dH*1000/Tb, Tb, dCp);

% MMH-H+
dS = dH*1000/Tb-30;
dH = 102.03; %kJ/mol

MakeLiquid('CH3NHNH3+',db.CH3NHNH2, dH, dS, 0, 0);

% CH3NHNH
dH = 37.21; %kJ/mol, Schmidt p 296
Tb = 365; %K
dCp = -45; %J/mol/K (based on hydrazine)

MakeLiquid('CH3NHNHL',db.CH3NHNH, dH, dH*1000/Tb, Tb, dCp);

% CH3NN
dH = 37.21; %kJ/mol, Schmidt p 296
Tb = 365; %K
dCp = -45; %J/mol/K (based on hydrazine)

MakeLiquid('CH3NN',db.CH3NN, dH, dH*1000/Tb, Tb, dCp);

% CH3NN+
dH = 102.03; %kJ/mol, Schmidt p 296
Tb = 365; %K
dS = dH*1000/Tb-30;

MakeLiquid('CH3NN+',db.CH3NN, dH, dS, 0, 0);


% CH3NHNH+
dH = 102.03; %kJ/mol, Schmidt p 296
Tb = 365; %K
dS = dH*1000/Tb-30;

MakeLiquid('CH3NHNH+',db.CH3NHNH, dH, dS, 0, 0);

% CH3N3
dH = 37.21; %kJ/mol, Schmidt p 296
Tb = 365; %K
dCp = -45; %J/mol/K (based on hydrazine)

MakeLiquid('CH3N3',db.CH3N3, dH, dH*1000/Tb, Tb, dCp);


%H2O
dH = 40.68; %kJ/mol, Schmidt p 296
Tb = 373.15; %K
dCp = -41.76; %J/mol/K

MakeLiquid('H2OL', db.H2O, dH, dH*1000/Tb, Tb, dCp);

%NO2
dH = 33.71; %kJ/mol, Schmidt p 296
Tb = 293.8; %K
dCp = 0; %J/mol/K (guessed)

MakeLiquid('NO2L',db.NO2, dH, dH*1000/Tb, Tb, dCp);

%HNO3
dH = 33.6; %kJ/mol
Tb = 330; %K
dCp = 0; %J/mol/K (guessed)

MakeLiquid('HNO3L',db.HNO3, dH, dH*1000/Tb, Tb, dCp);


%HONO
dH = 33.6; %kJ/mol, Schmidt p 296
Tb = 330; %K
dCp = 0; %J/mol/K (guessed)

MakeLiquid('HONOL',db.HONO, dH, dH*1000/Tb, Tb, dCp);


%CH3OH (methanol)
dH = 38.278; %kJ/mol, Schmidt p 296
Tb = 338; %K
dCp = -37; %J/mol/K (Burcat)

MakeLiquid('CH3OHL',db.CH3OH, dH, dH*1000/Tb, Tb, dCp);


%CH3ONO2 (methyl nitrate)
dH = 34.1; %NIST webbook
Tb = 338; %K, NIST webbook
dCp = 0; %J/mol/K

MakeLiquid('CH3ONO2L',db.CH3ONO2, dH, dH*1000/Tb, Tb, dCp);
