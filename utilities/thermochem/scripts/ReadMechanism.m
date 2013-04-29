clear all
close all
clc

% Load and compare reaction mechanisms from Chemkin files

Ns = 41;   %select 25, 29, or 41

chemfile = strcat('chemkin',filesep,'chem',num2str(Ns,'%2d'),'.inp');
thermofile = strcat('chemkin',filesep,'therm',num2str(Ns,'%2d'),'.dat');

T = 200:5:2500;
reactions = ReadCHEMKINReactions(chemfile);
species = ReadCHEMKINThermo(thermofile,T);

Ns =  length(fieldnames(species));
Nr = length(reactions);

specie_names = fieldnames(species);

for s = 1:Ns
    name = specie_names{s};
    is_used = false;
    
    for r = 1:Nr
        rxn = reactions(r);
        
        pn = {rxn.products.name};
        rn = {rxn.reactants.name};
        
        if any(strcmp(name, pn)) || any(strcmp(name, rn))
            is_used = true;
            break;
        end
    end
    
    if ~is_used
        fprintf('Specie %s is not used\n',name);
    end
end