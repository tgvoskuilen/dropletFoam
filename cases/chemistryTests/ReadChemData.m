clear all
close all
clc

% Variables:
%  p (1 bar)
%  Ti (280 K, 400 K)
%  N (1, 2, 4)

% N2O4 tried, but requires Ti over 800 K to initiate any reaction

files = dir('data/*.out');



data(size(files)) = struct();

for r = 1:length(files)
    filename = strcat('data',filesep,files(r).name);
    
    rawdata = importdata(filename,'\t',1);
    
    data(r).time = rawdata.data(:,1)*1e3; %ms
    data(r).T = rawdata.data(:,2); %K
    data(r).P = rawdata.data(1,3)/1e5; %bar
    
    C = textscan(filename(6:end),'%s T%f OF%f p%f constV %d.out','Delimiter','_');
    data(r).ox = C{1}{1};
    data(r).OF = C{3};
    data(r).Ns = C{5};
    data(r).Ti = C{2};
end


figure;

subplot(2,1,1)
hold all
set(gca,'Box','on','LineWidth',1,'FontSize',14)
xlabel('Time (ms)')
ylabel('T (K)')


for i = 1:length(data)
    if strcmp(data(i).ox,'NTO') && data(i).OF == 2
        if data(i).Ti == 300
            c = 'b';
        elseif data(i).Ti == 400
            c = 'k';
        elseif data(i).Ti == 600
            c = 'm';
        else
            c = 'r';
        end
        if data(i).Ns == 25
            l = '--';
        elseif data(i).Ns == 29
            l = ':';
        else
            l = '-';
        end
        plot(data(i).time, data(i).T,'Color',c,'LineStyle',l,'LineWidth',2);
    end
end
title('NTO');
xlim([0 0.5])

subplot(2,1,2)
hold all
set(gca,'Box','on','LineWidth',1,'FontSize',14)
xlabel('Time (ms)')
ylabel('T (K)')
title('RFNA')
for i = 1:length(data)
    if strcmp(data(i).ox,'RFNA') && data(i).OF == 2
        if data(i).Ti == 300
            c = 'b';
        elseif data(i).Ti == 400
            c = 'k';
        elseif data(i).Ti == 600
            c = 'm';
        else
            c = 'r';
        end
        if data(i).Ns == 25
            l = '--';
        elseif data(i).Ns == 29
            l = ':';
        else
            l = '-';
        end
        plot(data(i).time, data(i).T,'Color',c,'LineStyle',l,'LineWidth',2);
    end
end
xlim([0 0.5])

h300 = plot(-1,0,'-b','LineWidth',2);
h400 = plot(-1,0,'-k','LineWidth',2);
h600 = plot(-1,0,'-m','LineWidth',2);
h800 = plot(-1,0,'-r','LineWidth',2);
h25 = plot(-1,0,'--k','LineWidth',2);
h29 = plot(-1,0,':k','LineWidth',2);
h41 = plot(-1,0,'-k','LineWidth',2);

lh = legend([h300 h400 h600 h800 h25 h29 h41],...
    'T_i = 300 K',...
    'T_i = 400 K',...
    'T_i = 600 K',...
    'T_i = 800 K',...
    '25 specie',...
    '29 specie',...
    '41 specie');
set(lh,'Location','SouthOutside');

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'-dtiff','-r600','OxComparison.tif');