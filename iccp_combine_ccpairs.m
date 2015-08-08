function ccpairs_total = iccp_combine_ccpairs
% iccp_combine_ccpairs Combine all ccpairs data into struct array
%
% ccpairs_total = iccp_combine_ccpairs
% -----------------------------------------------------------
% Looks in the current folder for *-strfcmb-pairs-ccpairs.mat files, and
% then combines all the data in the files into the struct array
% ccpairs_total.


narginchk(0,0);

pfiles = dir('*-strfcmb-pairs-ccpairs.mat');
ccpairs_total = [];

for j = 1:length(pfiles)
   ccfile = pfiles(j).name;
   load(ccfile,'ccpairs');
   ccpairs_total = [ccpairs_total ccpairs];
   clear ccpairs;
end % (for j)

return;

