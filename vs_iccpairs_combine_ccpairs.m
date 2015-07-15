function ccpairs_total = vs_iccpairs_combine_ccpairs
% vs_iccpairs_spktype_strf_crosscorr Correlation between pairs
%
% vs_iccpairs_spktype_strf_crosscorr(override)
% -----------------------------------------------------------
%

%
% If the elements match, the cross-correlation between the spike trains for
% each pair are estimated, and the results are saved in a struct array that
% contains all the data in the individual files as well as the new
% cross-correlation data.
%
% The new struct array is then saved in a file with the ending: 
% *-strfcmb-pairs-ccpairs.mat
%
% If the override input argument is not supplied, then the function first
% checks to make sure that there is not already an output *-ccpairs.mat
% file. If there is, then the function skips processing the data for that
% set of files. If override = 1, then the function runs and writes over any
% *-ccpairs.mat file that is present.
%
% caa 8/15/14


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

