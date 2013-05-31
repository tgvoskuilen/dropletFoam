clear all
close all
clc


% Load thermodynamic database
db = ReadDb();
hz = db.N2H4;
dHf = H(hz,298);