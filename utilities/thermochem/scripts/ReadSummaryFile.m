clear all
close all
clc

filenames = {'summaryNoNADecomp41spMech.csv',...
             'summaryNADecomp41spMech.csv',...
             'summaryNoNADecompLiqRxn.csv',...
             'summaryNADecompLongerRun.csv'...
             'summaryReduced1.csv',...
             'summaryReduced2.csv'};

data = cell(length(filenames),1);
for i = 1:length(filenames)
    data{i} = importdata(filenames{i},',',1);

    data{i}.Ns = length(data{i}.textdata)-1;
    data{i}.t = data{i}.data(:,1).*1e3;



    figure;
    hold all
    h = [];
    for j = 1:data{i}.Ns
        if max(data{i}.data(:,j+1) < 1) && max(data{i}.data(:,j+1)) > 1e-6
            hN = plot(data{i}.t,data{i}.data(:,1+j));
            h = [h; hN]; %#ok<AGROW>
            set(hN,'DisplayName',data{i}.colheaders{j+1});
        end
    end
    xlabel('t (ms)')
    ylabel('Y')
    lh = legend(h,'Location','BestOutside');

    data{i}.Tidx = find(strcmpi(data{i}.colheaders,'T'));
    data{i}.T = data{i}.data(:,data{i}.Tidx);

end

figure;
hold all
for i = 1:length(data)
    plot(data{i}.t,data{i}.T)
end