function [id,table] = iccp_fracmb_to_table
% iccp_fracmb_to_table Tabulate all neuron FRA data parameters in folder
% 
%     [id,table] = iccp_fracmb_to_table
%
%     Goes through all '*-18115-fracmb-pairs-isi-fixed.mat' files in the
%     current folder, loads the fra struct array in each file. Each element
%     of fra holds the data for one neuron. The identifiers for each neuron
%     are added to a table and to a cell array, which can later be used to
%     determine how many neurons were recorded from.
% 
% 
%     table : struct of tables holding the data. Each row of one
%     of the tables holds the data for one neuron. Since only FRA data are
%     being processed, there is only one data table in table.
% 
%     id : cell array of string that identifies each neuron. The string 
%     includes experiment, site, channel, model, position, stimulus, atten, 
%     and nreps values.
% 
%     id may later be used to identify the unique neurons contained in 
%     ccpairs. To do this use the unique.m command in matlab.
% 
%         [c,ia,ic] = unique(id);
% 
%     and c(ia) will give the unique identifiers.
% 
%     Example: To get the unique values in a table:
% 
%         mtfUnique = table.mtf(ia,:);


d = dir('2004-*-fs18115-fracmb-pairs-isi-fixed.mat');
tableFRA = [];

id = {};

for ii = 1:length(d)

    filename = d(ii).name;
    load(filename, 'fra');
    fprintf('Processing %s\n', filename);

    for i = 1:length(fra)

        % Recording specifications
        exp = fra(i).exp;
        dn = datenum(exp);
        site = fra(i).site;
        chan = fra(i).chan;
        model1 = fra(i).model(1);
        position = fra(i).position;
        stim = fra(i).stim;
        index = findstr(stim, 'tc');
        stimNum = str2double(stim(index+2:end));
        atten = fra(i).atten;
        nreps = fra(i).nreps;


        % Put unit details in table
        tableFRA = [tableFRA; dn site chan model1 position stimNum atten nreps];


        % String identifier for each neuron in each pair
        id1 = sprintf('%s-%.0f-%.0f-%.0f-%s-%.0f-%.0f', ...
            exp, site, chan, model1, stim, atten, nreps);

        % Save identifiers to cell array
        id{length(id)+1} = id1;

    end % (for i)

end % (for ii)

FRACols = {'dn', 'site', 'chan', 'model', 'position', 'stimNum', ...
          'atten', 'nreps'};


% Save names and tables for each type of data
table.fraNames = FRACols;
table.fra = tableFRA;


return;







