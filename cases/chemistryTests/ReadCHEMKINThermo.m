function [species, T] = ReadCHEMKINThermo(filename,T)
% Reads a CHEMKIN formatted thermo file and calculates properties

%-------------------------------------------------------------------------
% Copyright (c) 2012, Tyler Voskuilen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR  
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
% CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%-------------------------------------------------------------------------
    
    fin = fopen(filename,'r');
    A = fscanf(fin,'%c');
    fclose(fin);

    line_ends = find(A==char(10));
    num_lines = length(line_ends);
    line_starts = [0 line_ends(1:end-1)];
    
    
    ns = 0;
    species(1) = struct();
    
    for l = 3:num_lines
        
        line = A(line_starts(l)+1:line_ends(l)-1);

        if length(line) == 80
            if isnan(str2double(line(1:15)))
                ns = ns+1;
                name = textscan(line(1:18),'%s');
                name = name{1}{1};
                
                % Clean parenthesis for solids (e.g. 'CH2(S)' -> 'CH2_S')
                name(name=='(') = '_';
                name(name==')') = '';
                
                species.(name).Tlow  = str2double(line(46:55));
                species.(name).Tmed  = str2double(line(66:73));
                species.(name).Thigh = str2double(line(56:65));
                
                species.(name).name = name; %#ok<*AGROW,*SAGROW>
                species.(name).a = [];
            else
                n(1) = str2double(line(1:15));
                n(2) = str2double(line(16:30));
                n(3) = str2double(line(31:45));
                n(4) = str2double(line(46:60));
                n(5) = str2double(line(61:75));
                
                if isnan(n(5))
                    n=n(1:4);
                end
                
                species.(name).a = [species.(name).a n];
            end
        end
    end
    
    spn = fieldnames(species);
    
    for s = 1:length(spn)
        species.(spn{s}).al = species.(spn{s}).a(8:14);
        species.(spn{s}).ah = species.(spn{s}).a(1:7);
    end

   
        
    
    
    dT = 10;
    R = 8.314; %J/mol-K
    
    if ~exist('T','var')
        T = 200:dT:5000;
    end

    for s = 1:length(spn)
        species.(spn{s}).Cp = Cp(species.(spn{s}),T) .* R;
        species.(spn{s}).h = H(species.(spn{s}),T) .* R;
        species.(spn{s}).s = S(species.(spn{s}),T) .* R;
        species.(spn{s}).g = species.(spn{s}).h - T.*species.(spn{s}).s;
        species.(spn{s}).hc = interp1(T, species.(spn{s}).h, 298);
        species.(spn{s}).hs = species.(spn{s}).h - species.(spn{s}).hc;
    end
    
end

function Cp = Cp(specie,T)

    % each specie has two sets of 7 coefficients
    lc = specie.al;
    hc = specie.ah;
    Tm = specie.Tmed;
    
    plc = lc(5:-1:1);
    phc = hc(5:-1:1);
    
    Cp = zeros(size(T));
    Cp(T<=Tm) = polyval(plc,T(T<=Tm));
    Cp(T>Tm) = polyval(phc,T(T>Tm));
end

function h = H(specie,T)

    % each specie has two sets of 7 coefficients
    lc = specie.al;
    hc = specie.ah;
    Tm = specie.Tmed;
    
    plc = [lc(5:-1:1)./(5:-1:1), lc(6)];
    phc = [hc(5:-1:1)./(5:-1:1), hc(6)];
    
    h = zeros(size(T));
    h(T<=Tm) = polyval(plc,T(T<=Tm));
    h(T>Tm) = polyval(phc,T(T>Tm));
end

function s = S(specie,T)

    % each specie has two sets of 7 coefficients
    lc = specie.al;
    hc = specie.ah;
    Tm = specie.Tmed;
    
    plc = [lc(5:-1:2)./(4:-1:1) lc(7)];
    phc = [hc(5:-1:2)./(4:-1:1) hc(7)];
    
    s = zeros(size(T));
    s(T<=Tm) = (lc(1)*log(T(T<=Tm)) + polyval(plc,T(T<=Tm)));
    s(T>Tm) = (hc(1)*log(T(T>Tm)) + polyval(phc,T(T>Tm)));
end
