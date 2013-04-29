clear all
close all
clc


% read LJ_Coeffs.csv and calculate sutherland coefficients As and Ts

filename = 'LJ_Coeffs.csv';

f = fopen(filename,'r');
C = textscan(f,'%s','Delimiter',',');
fclose(f);

ncols = 5;
N = (length(C{1})-ncols)/ncols;

data(1:N) = struct();
nameIDs = ncols+1:ncols:length(C{1});
sigmaIDs = nameIDs + 1;
TIDs = sigmaIDs + 1;
WIDs = TIDs + 1;
refIDs = WIDs + 1;

for i = 1:N
    data(i).name = C{1}{nameIDs(i)};
    data(i).ref = C{1}{refIDs(i)};
    data(i).sigma = str2double(C{1}{sigmaIDs(i)});
    data(i).T = str2double(C{1}{TIDs(i)});
    data(i).W = str2double(C{1}{WIDs(i)});
end

%Calculate Lennard Jones viscosity
T = 280:5:2000; %T range to fit to
for i = 1:N
    Ts = T ./ data(i).T;
    Omega = 1.16145.*Ts.^(-0.14874) + 0.52487.*exp(-0.77320.*Ts) + ...
        2.16178.*exp(-2.43787.*Ts);
    data(i).muLJ = 26.69e-5*sqrt(data(i).W.*T)./(data(i).sigma^2.*Omega); %Pa-s
end


%Load NASA database and compare
nasafile = 'trans_new2.txt';
f = fopen(nasafile,'r');
A = fscanf(f,'%c');
fclose(f);

line_ends = find(A==char(10));
num_lines = length(line_ends);
line_starts = [0 line_ends(1:end-1)];
species = struct();
ns = 0;

for l = 10:num_lines
    line = A(line_starts(l)+1:line_ends(l)-1);
    if ~isempty(line)
        if length(line) > 5
            start_char = line(1);
            if ~strcmp(start_char,' ')
                ns = ns + 1;
                name = textscan(line(1:10),'%s');
                name = name{1}{1};
                name(name=='(' | name==')' | name==',' | name=='-' | name=='!') = '_';
                if ~isfield(species,name)
                    species.(name).name = name;
                    species.(name).V = [];
                end
            else
                %species.(name).lines = [species.(name).lines; line];
                mode = line(2);
                Tl = str2double(line(3:10));
                Th = str2double(line(11:19));
                a(1) = str2double(line(22:35));
                a(2) = str2double(line(36:50));
                a(3) = str2double(line(51:65));
                a(4) = str2double(line(66:80));

                
                if strcmp(mode,'V')
                    species.(name).V = [species.(name).V; [Tl Th a]];
                end                
            end
        end
    end
end


for i = 1:N
    if isfield(species,data(i).name)
        s = size(species.(data(i).name).V);
        nTranges = s(1);
        
        A = species.(data(i).name).V(1,3);
        B = species.(data(i).name).V(1,4);
        C = species.(data(i).name).V(1,5);
        D = species.(data(i).name).V(1,6);
        data(i).muNASA1 = 1e-5.*exp(A.*log(T) + B./T + C./T.^2 + D);
        
    end
end


%plot comparison of muLJ and muNASA1
% for i = 1:N
%     if data(i).muNASA1
%         figure;
%         hold all
%         plot(T,data(i).muLJ,'b')
%         plot(T,data(i).muNASA1,'r');
%     end
% end


% Determine sutherland transport coefficients using muLJ
muSutherland = @(C,T) C(1).*T.^(1.5)./(T+C(2));
errMu = @(C,T,mu) sum(abs(muSutherland(C,T)-mu));

for i = 1:N
    mu0est = data(i).muLJ(1);
    data(i).sCoeffs = fminsearch(@(C) errMu(C,T,data(i).muLJ), [mu0est 300]);
    %data(i).muS = data(i).sCoeffs(1).*T.^1.5./(T+data(i).sCoeffs(2));
end


% Write Sutherland database output
fout = fopen('SutherlandCoeffs.csv','w');
fprintf(fout,'Specie,As,Ts\n');
for i = 1:N
    fprintf(fout,'%s,%7.6e,%7.6e\n',data(i).name,data(i).sCoeffs);
end
fclose(fout);

