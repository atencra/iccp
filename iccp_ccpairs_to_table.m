function [id,table] = iccp_ccpairs_to_table(ccpairs)
% iccp_ccpairs_to_table Tabulate all paired neuron data parameters
% 
%     [id,table] = iccp_ccpairs_to_table(ccpairs)
% 
%     ccpairs : struct array holding paired neuron data. Paired neuron data is 
%     data for two neurons recorde from the same electrode channel.
% 
%     The data normally includes STRF firing rate, phase-locking index, 
%     cf, q, latency, bmf, and  nonlinearity parameters.
% 
%     table : struct of tables holding the paired data. Each row of one
%     of the tables holds the data for one member of a pair of neurons.
%     Since each element of ccpairs holds data for two neurons, each table
%     will have 2*length(ccpairs) rows.
% 
%     id : cell array of string that identifies each neuron in a pair in
%     ccpairs. The string includes experiment, site, channel, model,
%     position, stimulus, and atten value. length(id) = 2*length(ccpairs).
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




tablePT = [];
tableParams = [];
tableMTF = [];
tableFIO = [];

id = {};

for i = 1:length(ccpairs)

    % Recording specifications
    exp = ccpairs(i).exp;
    dn = datenum(exp);
    site = ccpairs(i).site;
    chan = ccpairs(i).chan;
    model1 = ccpairs(i).model1(1);
    model2 = ccpairs(i).model2(1);
    position = ccpairs(i).position;
    stim = ccpairs(i).stim;
    index = findstr(stim, 'dmr');
    stimNum = str2double(stim(index+3:end));
    atten = ccpairs(i).atten;


    % Pure tone params
    cf1 = ccpairs(i).cf1;
    cf2 = ccpairs(i).cf2;
    q1 = ccpairs(i).q1;
    q2 = ccpairs(i).q2;
    lat1 = ccpairs(i).latency1;
    lat2 = ccpairs(i).latency2;
    tablePT = [tablePT; dn site chan model1 position stimNum atten cf1 q1 lat1];
    tablePT = [tablePT; dn site chan model2 position stimNum atten cf2 q2 lat2];
    PTCols = {'dn', 'site', 'chan', 'model', 'position', 'stimNum', ...
              'atten', 'cf', 'q', 'lat'};


    % STRF FR and PLI
    fr1 = ccpairs(i).fr1;
    fr2 = ccpairs(i).fr2;
    pli1 = ccpairs(i).pli1;
    pli2 = ccpairs(i).pli2;
    tableParams = [tableParams; dn site chan model1 position stimNum atten fr1 pli1];
    tableParams = [tableParams; dn site chan model1 position stimNum atten fr2 pli2];
    ParamsCols = {'dn', 'site', 'chan', 'model', 'position', 'stimNum', ...
              'atten', 'fr', 'pli'};


    % STRF BMFs
    btmf1 = ccpairs(i).btmf1;
    btmf2 = ccpairs(i).btmf2;
    bsmf1 = ccpairs(i).bsmf1;
    bsmf2 = ccpairs(i).bsmf2;
    tableMTF = [tableMTF; dn site chan model1 position stimNum atten btmf1 bsmf1];
    tableMTF = [tableMTF; dn site chan model1 position stimNum atten btmf2 bsmf2];
    MTFCols = {'dn', 'site', 'chan', 'model', 'position', 'stimNum', ...
              'atten', 'btmf', 'bsmf'};


    % STRF BMFs
    fioasi1 = ccpairs(i).fioasi1;
    theta1 = ccpairs(i).theta1;
    sigma1 = ccpairs(i).sigma1;
    nmse1 = ccpairs(i).nmse1;
    r21 = ccpairs(i).r21;

    fioasi2 = ccpairs(i).fioasi2;
    theta2 = ccpairs(i).theta2;
    sigma2 = ccpairs(i).sigma2;
    nmse2 = ccpairs(i).nmse2;
    r22 = ccpairs(i).r22;

    tableFIO = [tableFIO; dn site chan model1 position stimNum atten ...
                    fioasi1 theta1 sigma1 nmse1 r21];

    tableFIO = [tableFIO; dn site chan model1 position stimNum atten ...
                    fioasi2 theta2 sigma2 nmse2 r22];

    FIOCols = {'dn', 'site', 'chan', 'model', 'position', 'stimNum', ...
              'atten', 'fioasi', 'theta', 'sigma', 'nmse', 'r2'};



    % String identifier for each neuron in each pair
    id1 = sprintf('%s-%.0f-%.0f-%.0f-%s-%.0f', ...
        exp, site, chan, model1, stim, atten);

    id2 = sprintf('%s-%.0f-%.0f-%.0f-%s-%.0f', ...
        exp, site, chan, model2, stim, atten);


    % Save identifiers to struct arrays
    id{2*(i-1)+1} = id1;
    id{2*i} = id2;

end % (for i)

% Save names and tables for each type of data
table.ptNames = PTCols;
table.pt = tablePT;
table.paramsNames = ParamsCols;
table.params = tableParams;
table.mtfNames = MTFCols;
table.mtf = tableMTF;
table.fioNames = FIOCols;
table.fio = tableFIO;



return;













