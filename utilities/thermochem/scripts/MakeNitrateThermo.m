clear all
close all
clc

% Load the thermo db
db = ReadDb();

name = 'MMH3N';
Ln = db.MMHN.low + 2*db.HNO3L.low;
Hn = db.MMHN.high + 2*db.HNO3L.high;

% Ln(6) = Ln(6) - dH/R;
% Hn(6) = Hn(6) - dH/R;
% 
% % dG = dH - T*dS
% dS = (dH - dG)*1000/298; %cal/mol/K
% 
% Ln(7) = Ln(7) - dS/R;
% Hn(7) = Hn(7) - dS/R;

fprintf('CHEMKIN Polynomial Coefficients for %s\n',name);
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    2\n',Hn(1),Hn(2),Hn(3),Hn(4),Hn(5));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E%+10.8E    3\n',Hn(6),Hn(7),Ln(1),Ln(2),Ln(3));
fprintf('%+10.8E%+10.8E%+10.8E%+10.8E                   4\n\n',Ln(4),Ln(5),Ln(6),Ln(7));