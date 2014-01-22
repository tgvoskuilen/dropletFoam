clear all
close all
clc

% Calculate equilibrium mass fractions of a mixture of NO2/N2O4 at a given
% temperature in different oxidizers


% Inputs, temperature and oxidizer
Ti=300;
ox = 'NTO';



T = Ti-20:1:Ti+20;
species = ReadCHEMKINThermo('oxChem.dat',T);

dH = species.('N2O4').h - 2*species.('NO2').h;
dS = species.('N2O4').s - 2*species.('NO2').s;

dG = dH - T.*dS;
Kp = exp(-dG./(8.314.*T));
calc_eq = true;

if strcmpi(ox,'nto')
    %NTO
    Y_HNO3 = 0;
    Y_NO2_orig = 1;
    Y_H2O = 0;
elseif strcmpi(ox,'rfna')
    % RFNA
    Y_HNO3 = 0.813;
    Y_NO2_orig = 0.187;
    Y_H2O = 0;
elseif strcmpi(ox,'wfna')
    % WFNA
    Y_HNO3 = 0.978;
    Y_NO2_orig = 0.022;
    Y_H2O = 0;
elseif strcmpi(ox,'no2')
    Y_HNO3 = 0;
    Y_NO2_orig = 1;
    Y_H2O = 0;
    calc_eq = false;
else
    error('Invalid oxidizer selection')
end

% kgi/kmoli
W_HNO3 = 63;
W_NO2 = 46;
W_N2O4 = 92;
W_H2O = 18;
W_MMH = 46;


% kmoli/kgT
N_HNO3 = Y_HNO3/W_HNO3;
N_NO2_orig = Y_NO2_orig/W_NO2;

% kmoli/kmolT
x_HNO3 = @(N_N2O4) (N_HNO3 ./ (N_HNO3 + N_NO2_orig - N_N2O4));
x_N2O4 = @(N_N2O4) (N_N2O4 ./ (N_HNO3 + N_NO2_orig - N_N2O4));
x_NO2 = @(N_N2O4)  ((N_NO2_orig - 2.*N_N2O4) ./ (N_HNO3 + N_NO2_orig - N_N2O4));

err = @(N_N2O4) abs(Kp(T==Ti) .* x_NO2(N_N2O4).^2 - x_N2O4(N_N2O4));

if calc_eq
    n = fminbnd(err, 0, N_NO2_orig/2);
else
    n = 0;
end

% oxidizer molar mass
W_ox = x_HNO3(n)*W_HNO3 + x_NO2(n)*W_NO2 + x_N2O4(n)*W_N2O4;

% Calculate mass fractions
fprintf('Y_NO2 = %6.5f\n', x_NO2(n)*W_NO2/W_ox);
fprintf('Y_N2O4 = %6.5f\n', x_N2O4(n)*W_N2O4/W_ox);
fprintf('Y_HNO3 = %6.5f\n', x_HNO3(n)*W_HNO3/W_ox);



