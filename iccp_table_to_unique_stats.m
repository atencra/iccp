function [tablePairsUnq] = iccp_table_to_unique_stats(id, table, nsamples)
% iccp_table_to_unique_stats Statistics for unique data table values
% 
%     tablePairsUnq = iccp_table_to_unique_stats(id, table, nsamples)
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
%
%     table : struct of tables holding the paired data. Each row of one
%     of the tables holds the data for one member of a pair of neurons.
%     Since each element of ccpairs holds data for two neurons, each table
%     will have 2*length(ccpairs) rows.
% 
%     The data normally includes STRF firing rate, phase-locking index, 
%     cf, q, latency, bmf, and  nonlinearity parameters.
%
%     ***** id and table are obtained from iccp_ccpairs_to_table.com
%
%     nsamples : if specified, samples the unique value nsamples # of times
%     and returns pairwise distributions for fr, pli, btmf, bsmf, cf, q,
%     latency, asi, theta, and sigma. nsamples should be set to
%     length(ccpairs) if it is used. The default value is 1029.
%
%     tablePairsUnq : struct with the randomly sampled pairs. Has fields of
%     firing rate, and other variables. Each field holds an nsamples X 2
%     array, where each row is one pair.


if nargin == 2
    nsamples = 1029;
end

[c, ia, ic] = unique(id);
length(c)
idUnq = id(ia);


% Get stats for strf fr and pli
indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

fr1 = table.params(indexSamp1,end-1);
fr2 = table.params(indexSamp2,end-1);
fr12 = [fr1(:) fr2(:)];


indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

pli1 = table.params(indexSamp1,end);
pli2 = table.params(indexSamp2,end);
pli12 = [pli1(:) pli2(:)];



% Get stats for strf pure tone params
indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

q1 = table.pt(indexSamp1,end-1);
q2 = table.pt(indexSamp2,end-1);
q12 = [q1(:) q2(:)];


indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

lat1 = table.pt(indexSamp1,end);
lat2 = table.pt(indexSamp2,end);
lat12 = [lat1(:) lat2(:)];



% Get stats for strf bmf

indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

btmf1 = table.mtf(indexSamp1,end-1);
btmf2 = table.mtf(indexSamp2,end-1);
btmf12 = [btmf1(:) btmf2(:)];

indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

bsmf1 = table.mtf(indexSamp1,end);
bsmf2 = table.mtf(indexSamp2,end);
bsmf12 = [bsmf1(:) bsmf2(:)];



% Get stats for fio params: asymmetry index

indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

fioasi1 = table.fio(indexSamp1,end-4);
fioasi2 = table.fio(indexSamp2,end-4);
fioasi12 = [fioasi1(:) fioasi2(:)];



% Get stats for fio params: theta, sigma. Only resample once, since
% the theta/sigma pair is linked.

indexSamp1 = ceil( nsamples * rand(1,nsamples) );
indexSamp2 = ceil( nsamples * rand(1,nsamples) );

theta1 = table.fio(indexSamp1,end-3);
theta2 = table.fio(indexSamp2,end-3);
theta12 = [theta1(:) theta2(:)];

sigma1 = table.fio(indexSamp1,end-2);
sigma2 = table.fio(indexSamp2,end-2);
sigma12 = [sigma1(:) sigma2(:)];

nmse1 = table.fio(indexSamp1,end-1);
nmse2 = table.fio(indexSamp2,end-1);
nmse12 = [nmse1(:) nmse2(:)];

r1 = table.fio(indexSamp1,end);
r2 = table.fio(indexSamp2,end);
r12 = [r1(:) r2(:)];


tablePairsUnq.fr12 = fr12;
tablePairsUnq.pli12 = pli12;
tablePairsUnq.q12 = q12;
tablePairsUnq.lat12 = lat12;
tablePairsUnq.btmf12 = btmf12;
tablePairsUnq.bsmf12 = bsmf12;
tablePairsUnq.fioasi12 = fioasi12;
tablePairsUnq.theta12 = theta12;
tablePairsUnq.sigma12 = sigma12;
tablePairsUnq.nmse12 = nmse12;
tablePairsUnq.r12 = r12;


% tablePT = [];
% tableParams = [];
% tableMTF = [];
% tableFIO = [];
% 
% tablePT = [tablePT; dn site chan model1 position stimNum atten cf1 q1 lat1];
% tableParams = [tableParams; dn site chan model1 position stimNum atten fr1 pli1];
% tableMTF = [tableMTF; dn site chan model1 position stimNum atten btmf1 bsmf1];
% tableFIO = [tableFIO; dn site chan model1 position stimNum atten ...
%                 fioasi1 theta1 sigma1 nmse1 r21];
% 
% table.pt = tablePT;
% table.params = tableParams;
% table.mtf = tableMTF;
% table.fio = tableFIO;



return;













