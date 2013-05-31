clear all
close all
clc

%TODO: Use L(T) expression
% L(T) = Lb*((Tc-T)/(Tc-Tb))^0.38 to get Cp change

% Load thermodynamic database
db = ReadDb();
R = 1.9858775; %cal/mol/K

% MMH Nitrite salt
dH = -51.5; %kcal/mol
dG = -2.1; %kcal/mol

T = 300:5:3000;

name = 'CH3NHNH3-ONO';
Ln = db.CH3NHNH2.low + db.HONO.low;
Hn = db.CH3NHNH2.high + db.HONO.high;

Ln(6) = Ln(6) - dH/R;
Hn(6) = Hn(6) - dH/R;

% dG = dH - T*dS
dS = (dH - dG)*1000/298; %cal/mol/K

Ln(7) = Ln(7) - dS/R;
Hn(7) = Hn(7) - dS/R;

fprintf('CHEMKIN Polynomial Coefficients for %s\n',name);
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    2\n',Hn(1),Hn(2),Hn(3),Hn(4),Hn(5));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    3\n',Hn(6),Hn(7),Ln(1),Ln(2),Ln(3));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E                   4\n\n',Ln(4),Ln(5),Ln(6),Ln(7));