clear all
close all
clc

files = dir('data2/*.csv');
data(1:length(files)) = struct();

for i = 1:length(files)
    raw = importdata(fullfile('data2',files(i).name),',',1);
    C = textscan(files(i).name,'%4s_T%f_OF%f_p%f_constV_%f.csv');
    
    data(i).ox = C{1}{1};
    data(i).Ti = C{2};
    data(i).OF = C{3};
    data(i).Pi = C{4};
    data(i).Nspecies = C{5};
    
    %Get arrays from raw
    for j = 1:length(raw.colheaders)
        data(i).(raw.colheaders{j}) = raw.data(:,j);
    end
    
    data(i).Tmax = max(data(i).T);
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

Tis = sort(unique([data.Ti]));
colors = hsv(length(Tis));

figure;
hold on
h = zeros(size(data));
for i = 1:length(data)
    h(i) = plot(data(i).OF, data(i).Tmax,'Line','None','LineWidth',1.5);
    
    % set marker based on oxidizer
    if strcmpi(data(i).ox,'rfna')
        set(h(i),'Marker','s');
    elseif strcmpi(data(i).ox,'wfna')
        set(h(i),'Marker','o');
    else
        %nothing
    end
    
    % set color based on initial temperature
    cid = find(Tis==data(i).Ti,1,'first');
    if data(i).Nspecies == 25
        set(h(i),'Color',colors(cid,:),'MarkerFaceColor',colors(cid,:));
    elseif data(i).Nspecies == 41
        set(h(i),'Color',colors(cid,:),'MarkerFaceColor','w');
    end
    
    
    
    
end




