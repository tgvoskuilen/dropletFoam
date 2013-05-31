function data = ReadDb()

    filename = '../props/thermo.pydb';

    fid = fopen(filename,'r');
    
    %list of molar masses for atoms used in here
    Ws.C = 12.011;
    Ws.N = 14.007;
    Ws.H = 1.008;
    Ws.O = 15.999;
    Ws.AR = 39.948;
    Ws.HE = 4.0026;
    Ws.E = 0;
    
    data = struct();

    line = fgets(fid);
    while ischar(line)
        
        A = textscan(line,'%s');
        entry = struct();
        for j = 1:length(A{1})
            ss = A{1}{j};

            if strcmp(ss,'''Thigh'':')
                entry.Thigh = str2double(A{1}{j+1}(1:end-1));
            end

            if strcmp(ss,'''Tlow'':')
                entry.Tlow = str2double(A{1}{j+1}(1:end-1));
            end

            if strcmp(ss,'''Tmid'':')
                entry.Tmid = str2double(A{1}{j+1}(1:end-1));
            end

            if strcmp(ss,'''Ts'':')
                entry.Ts = str2double(A{1}{j+1}(1:end-1));
            end

            if strcmp(ss,'''As'':')
                entry.As = str2double(A{1}{j+1}(1:end-1));
            end

            if ~isempty(strfind(ss,'lowCpCoeffs'))
                entry.low(1) = str2double(A{1}{j+1}(2:end-1));
                entry.low(2) = str2double(A{1}{j+2}(1:end-1));
                entry.low(3) = str2double(A{1}{j+3}(1:end-1));
                entry.low(4) = str2double(A{1}{j+4}(1:end-1));
                entry.low(5) = str2double(A{1}{j+5}(1:end-1));
                entry.low(6) = str2double(A{1}{j+6}(1:end-1));
                entry.low(7) = str2double(A{1}{j+7}(1:end-2));
            end

            if ~isempty(strfind(ss,'highCpCoeffs'))
                entry.high(1) = str2double(A{1}{j+1}(2:end-1));
                entry.high(2) = str2double(A{1}{j+2}(1:end-1));
                entry.high(3) = str2double(A{1}{j+3}(1:end-1));
                entry.high(4) = str2double(A{1}{j+4}(1:end-1));
                entry.high(5) = str2double(A{1}{j+5}(1:end-1));
                entry.high(6) = str2double(A{1}{j+6}(1:end-1));
                entry.high(7) = str2double(A{1}{j+7}(1:end-2));
            end

            if strcmp(ss,'''Name'':')
                entry.name = A{1}{j+1}(2:end-2);
            end

            if ~isempty(strfind(ss,'Atoms'))
                k = j+1;
                atoms = struct();
                while true
                    if strfind(A{1}{k},'{')
                        L = A{1}{k}(3:end-2);
                    else
                        L = A{1}{k}(2:end-2);
                    end


                    if strfind(A{1}{k+1},'}')
                        n = str2double(A{1}{k+1}(1:end-2));
                        kill = true;
                    else
                        n = str2double(A{1}{k+1}(1:end-1));
                        kill = false;
                    end

                    atoms.(L) = n;

                    if kill
                        break;
                    else
                        k = k + 2;
                    end
                end
                entry.atoms = atoms;
            end
        end
        
        f = fieldnames(entry.atoms);
        W = 0;

        for a = 1:length(f)
            W = W + Ws.(f{a})*entry.atoms.(f{a});
        end
        entry.W = W;

        v = genvarname(entry.name);
        data.(v) = entry;

        line = fgets(fid);
    end
    fclose(fid);
    