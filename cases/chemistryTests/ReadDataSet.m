clear all
close all
clc

folder = 'data';

files = dir([folder,filesep,'*.csv']);
data = struct();
i = 1;
for j = 1:length(files)
    raw = importdata(fullfile(folder,files(j).name),',',1);
    C = textscan(files(j).name,'%4s_T%f_OF%f_p%f_constP_%f.csv');

    if isempty(C{5})
        C = textscan(files(j).name,'%4s_T%f_OF%f_p%f_constP_v2_%f.csv');
        if isempty(C{5})
            continue
        end
        
        C{5} = 42;
    end
    
    data(i).ox = C{1}{1};
    data(i).Ti = C{2};
    data(i).OF = C{3};
    data(i).Pi = C{4};
    data(i).Nspecies = C{5};
    
    switch data(i).Nspecies
        case 25
            data(i).mechName = 'RChem1';
        case 29
            data(i).mechName = 'RChem2';
        case 41
            data(i).mechName = 'RChem3';
        case 42
            data(i).mechName = 'RChem3S';
        otherwise
            data(i).mechName = 'NaN';
    end

    %Get arrays from raw
    for k = 1:length(raw.colheaders)
        hn = raw.colheaders{k};
        hn(hn=='(' | hn==')') = '_';
        data(i).(hn) = raw.data(:,k);
    end

    data(i).Tmax = max(data(i).T);
    data(i).Tend = data(i).T(end);

    i = i + 1;
end

% make plot of Tmax vs O/F
% rfna.data = data(strcmpi({data.ox},'rfna'));
% wfna.data = data(strcmpi({data.ox},'wfna'));
% 
% wfna.OF = [wfna.data.OF];
% wfna.Tmax = [wfna.data.Tmax];
% 
% rfna.OF = [rfna.data.OF];
% rfna.Tmax = [rfna.data.Tmax];

load('CEA.mat');

Tis = sort(unique([data.Ti]));

%% Plot T vs OF
figure;
set(gcf,'Units','inches','Position',[2 7 3.5 2.75]);
hold on

rfnaColor = 'k';
wfnaColor = [.75 .75 .75];

plot(CEA.OF, CEA.T_rfna,'Line','-','Color',rfnaColor,'LineWidth',2)
plot(CEA.OF, CEA.T_wfna,'Line','-','Color',wfnaColor,'LineWidth',2)

subset = data([data.Nspecies] ~= 29);

h = zeros(size(subset));
for i = 1:length(subset)
    h(i) = plot(subset(i).OF, subset(i).Tmax,'Line','None',...
        'LineWidth',1,'MarkerSize',7);
    %plot(data(i).OF, data(i).Tend,'xk');
    
    % set marker based on oxidizer
    if strcmpi(subset(i).ox,'rfna')
        set(h(i),'Color','k','MarkerFaceColor',rfnaColor);
    elseif strcmpi(subset(i).ox,'wfna')
        set(h(i),'Color','k','MarkerFaceColor',wfnaColor);
    else
        %nothing
    end
    
    % set color based on initial temperature
    %cid = find(Tis==data(i).Ti,1,'first');
    if subset(i).Nspecies == 25
        set(h(i),'Marker','s');
    elseif subset(i).Nspecies == 41
        set(h(i),'Marker','o'); % full RChem3
    elseif subset(i).Nspecies == 42
        set(h(i),'Marker','d'); % stabilized RChem3
    end
    
end
xlabel('Mass O/F (Fuel = MMH)')
ylabel('T (K)')
set(gca,'Box','on','LineWidth',1)
xlim([0 8])
%hRC2 = plot(20,600,'^k','LineWidth',1,'MarkerSize',7,'DisplayName','RChem2');
hRC3 = plot(20,600,'ok','LineWidth',1,'MarkerSize',7,'DisplayName','RChem2');
hRC3s = plot(20,600,'dk','LineWidth',1,'MarkerSize',7,'DisplayName','RChem2S');

hRC1 = plot(20,600,'sk','LineWidth',1,'MarkerSize',7,'Displayname','RChem1');
hRFNA = plot(20,600,'sk','Color',rfnaColor,'MarkerFaceColor',rfnaColor,'MarkerSize',9,'DisplayName','RFNA');
hWFNA = plot(20,600,'sk','Color',wfnaColor,'MarkerFaceColor',wfnaColor,'MarkerSize',9,'DisplayName','WFNA');
hCEA = plot([20 21],[600 600],'-k','LineWidth',2,'DisplayName','CEA');
ax1 = gca;
text(6,2700,'P = 1 bar')
ax2 = axes('Position',get(ax1,'Position'),'Visible','off');

lh1 = legend(ax1,[hRC1 hRC3 hRC3s],'Location','SouthEast');
lh2 = legend(ax2,[hWFNA hRFNA hCEA],'Location','SouthEast');
lh1pos = get(lh1,'Position');
lh2pos = get(lh2,'Position');
set(lh2,'Position',[lh1pos(1)-lh2pos(3)+.02, lh1pos(2), lh2pos(3) lh1pos(4)]);
set([lh1 lh2],'Box','off');

set(gcf,'PaperPositionMode','auto');
print(gcf,'-dtiff','-r600','Chem_CEA_Comparison.tif');


%% Plot T vs OF - just RChem2
figure;
set(gcf,'Units','inches','Position',[2 7 4 3]);
hold on

data2 = data(strcmpi({data.ox},'rfna'));

mechColors = {'r','c','g','y'};
mechMarkers = {'s','^','o','d'};
hCEA = plot(CEA.OF, CEA.T_rfna,'Line','-','Color','k','LineWidth',2,'DisplayName','CEA');

chem{1} = data2([data2.Nspecies]==25);
chem{2} = data2([data2.Nspecies]==29);
chem{3} = data2([data2.Nspecies]==41);
chem{4} = data2([data2.Nspecies]==42);

h = zeros(size(chem));
for i = [1 3 4 2]
    h(i) = plot([chem{i}.OF],[chem{i}.Tmax],...
        'Line','--',...
        'Marker',mechMarkers{i},...
        'LineWidth',1,...
        'MarkerSize',7,...
        'Color','k',...
        'MarkerFaceColor',mechColors{i},...
        'DisplayName',chem{i}(1).mechName);
end


xlabel('Mass O/F (RFNA/MMH)')
ylabel('T (K)')
set(gca,'Box','on','LineWidth',1)
xlim([0 8])
% hRC1 = plot(20,600,'sk','LineWidth',1,'MarkerSize',7,'Displayname','RChem1','MarkerFaceColor',mechColors{1});
% hRC2 = plot(20,600,'^k','LineWidth',1,'MarkerSize',7,'DisplayName','RChem2','MarkerFaceColor',mechColors{2});
% hRC3 = plot(20,600,'ok','LineWidth',1,'MarkerSize',7,'DisplayName','RChem3','MarkerFaceColor',mechColors{3});
% hRC3s = plot(20,600,'dk','LineWidth',1,'MarkerSize',7,'DisplayName','RChem3S','MarkerFaceColor',mechColors{4});
% 
% hCEA = plot([20 21],[600 600],'-k','LineWidth',2,'DisplayName','CEA');
ax1 = gca;
text(6,2700,'P = 1 bar')

lh1 = legend(ax1,[hCEA h],'Location','SouthEast');
% lh1pos = get(lh1,'Position');
% lh2pos = get(lh2,'Position');
% set(lh2,'Position',[lh1pos(1)-lh2pos(3)+.02, lh1pos(2), lh2pos(3) lh1pos(4)]);
set(lh1,'Box','off');

set(gcf,'PaperPositionMode','auto');
print(gcf,'-dtiff','-r600','Chem_CEA_Comparison_RFNA_Only.tif');



%% Make plots of T vs time
figure;
set(gcf,'Units','inches','Position',[2 2 12 4]);
OFlist = [0.5, 1, 2, 3, 4, 5, 6];
colors = hsv(length(OFlist));

ax(1) = subplot(1,4,1);
hold all
ylim([500 3000])
ax(2) = subplot(1,4,2);
hold all
ylim([500 3000])
ax(3) = subplot(1,4,3);
hold all
ax(4) = subplot(1,4,4);
hold all
ylim([500 3000])

for i = 1:length(data)
    OFid = find(data(i).OF==OFlist,1,'first');
    
    if ~isempty(OFid)
        if data(i).Nspecies==25
            h=plot(ax(1),data(i).Time.*1000, data(i).T,'LineWidth',2);
        elseif data(i).Nspecies==29
            h=plot(ax(2),data(i).Time.*1000, data(i).T,'LineWidth',2);
        elseif data(i).Nspecies==41
            h=plot(ax(3),data(i).Time.*1000, data(i).T,'LineWidth',2);
        elseif data(i).Nspecies==42
            h=plot(ax(4),data(i).Time.*1000, data(i).T,'LineWidth',2);
        end
        
        if strcmpi(data(i).ox,'wfna')
            set(h,'Line','-');
        elseif strcmpi(data(i).ox,'rfna')
            set(h,'Line','--');
        end
        
        set(h,'Color',colors(OFid,:));
    end
end
title(ax(1),'RChem1');
title(ax(2),'RChem2');
title(ax(3),'RChem3');
title(ax(4),'RChem3S');
set(ax,'Box','on','LineWidth',1,'XScale','log');
for i = 1:4
    xlabel(ax(i),'Time (ms)')
    ylabel(ax(i),'T (K)')
end

%% Make plot of key species at O/F = 2

OFs = [data.OF];
oxs = {data.ox};
Nss = [data.Nspecies];

subset = data(OFs==2 & strcmpi(oxs,'rfna') & (Nss==25 | Nss==42 | Nss==29));


figure;
set(gcf,'Units','inches','Position',[1 1 12 3]);
ax(1) = subplot(1,3,1);
hold all
ax(2) = subplot(1,3,2);
hold all
ax(3) = subplot(1,3,3);
hold all
h = cell(3,1);
lh = zeros(3,1);

for s = 1:length(subset)
    title(ax(s), subset(s).mechName)
    flds = fieldnames(subset(s));
    
    noplot = {'Time','T','Rspecific','p','rho','phi'};
    yesplot = {'CH3NHNH2','CH3NNH','HNO3','CO','CO2','NO','HONO','N2','NO2','OH'};
    names = {'CH_3NHNH_2','CH_3NNH','HNO_3','CO','CO_2','NO','HONO','N_2','NO_2','OH'};
    h{s} = [];
    
    
    for n = 1:length(flds)
        if length(subset(s).(flds{n})) == length(subset(s).Time) && ...
                ~any(strcmpi(noplot, flds{n})) && ...
                any(strcmpi(yesplot, flds{n}))
            
            matchloc = find(strcmpi(yesplot, flds{n}),1,'first');
            
            hNew = plot(ax(s),subset(s).Time.*1000, subset(s).(flds{n}),'LineWidth',1.5);
            set(hNew,'DisplayName',names{matchloc});
            h{s} = [h{s}; hNew];
        end
    end
    set(h{s}(8:end),'Line','--');
    
    if s == 3
    lh(s) = legend(h{s},'Location','BestOutside');
    end
    xlabel(ax(s),'Time (ms)')
    ylabel(ax(s),'Mass Fraction');
    ylim(ax(s),[0 .6])
end



set(ax,'Box','on','LineWidth',1,'XScale','log');







