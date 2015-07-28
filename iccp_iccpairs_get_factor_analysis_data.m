function [mat,label] = vs_iccpairs_get_factor_analysis_data(ccpairs)
% vs_iccpairs_get_factor_analysis_data Get table for factor analysis
% 
%    [mat,label] = vs_iccpairs_get_factor_analysis_data(ccpairs) gets the data
%    from ccpairs for ccc,  cf,  q,  latency,  fr,  rpi, strfsi, fiosi, and asi
%    and calculates the difference between the values for each pair
%    of neurons.
% 
%    ccpairs is a struct array, where each element holds the data for one
%    pair of ICC neurons.
% 
%    The differences between parameter values are returned as columns 
%    of mat. Each row of mat is one pair of neurons.
%
%    label is a cell array of strings that denotes the datatype of each
%    column of mat.



ccc_pos = [ccpairs.ccc_pos];
pd_pos = [ccpairs.pd_pos];
hw_pos = [ccpairs.hw_pos];
sigpos = [ccpairs.sigfeature_pos];
index = find(sigpos < 1);
ccc_pos(index) = nan; % set nonsignificant values to nan

ccc_neg = [ccpairs.ccc_neg];
pd_neg = [ccpairs.pd_neg];
hw_neg = [ccpairs.hw_neg];
signeg = [ccpairs.sigfeature_neg];
index = find(signeg < 1);
ccc_neg(index) = nan; % set nonsignificant values to nan



cf1 = [ccpairs.cf1];
cf2 = [ccpairs.cf2];
cfdiff = abs( log2(cf1 ./ cf2) );

q1 = [ccpairs.q1];
q2 = [ccpairs.q2];
qdiff = abs( log2(q1 ./ q2) );

latency1 = [ccpairs.latency1];
latency2 = [ccpairs.latency2];
latencydiff = abs( latency1 - latency2);


fr1 = [ccpairs.fr1];
fr2 = [ccpairs.fr2];
frdiff = abs( fr1 - fr2);

rpi1 = [ccpairs.pli1];
rpi2 = [ccpairs.pli2];
rpidiff = abs( rpi1 - rpi2);

strfsi = [ccpairs.strfsimilarity];
strfsidiff = 1 - strfsi;

fiosi = [ccpairs.fiosimilarity];
fiosidiff = 1 - fiosi;


asi1 = [ccpairs.fioasi1];
asi2 = [ccpairs.fioasi2];
asidiff = abs( asi1 - asi2 );

theta1 = [ccpairs.theta1];
sigma1 = [ccpairs.sigma1];
nmse1 = [ccpairs.nmse1];

theta2 = [ccpairs.theta2];
sigma2 = [ccpairs.sigma2];
nmse2 = [ccpairs.nmse2];

thetadiff = abs(theta1 - theta2);
sigmadiff = abs(sigma1 - sigma2);
index = find(nmse1 > 0.2 | nmse2 > 0.2);
thetadiff(index) = nan;
sigmadiff(index) = nan;


% Can now add transition and threshold differences to this ...

mat = [ccc_pos(:) pd_pos(:) hw_pos(:) ccc_neg(:) pd_neg(:) hw_neg(:) ...
cfdiff(:) qdiff(:) latencydiff(:) frdiff(:) rpidiff(:) ...
strfsidiff(:) fiosidiff(:) asidiff(:) thetadiff(:) sigmadiff(:)];

label = {'ccc_pos',  'pd_pos', 'hw_pos', 'ccc_neg',  'pd_neg', 'hw_neg', ...
'cfdiff', 'qdiff', 'latencydiff', 'frdiff', 'rpidiff', ... 
'strfsidiff', 'fiosidiff', 'asidiff', 'thetadiff', 'sigmadiff'};


return;



