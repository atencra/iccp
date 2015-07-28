function fioFit_total = iccp_combine_fio_fitparams
% iccp_spktype_strf_crosscorr Correlation between pairs
%
% iccp_spktype_strf_crosscorr(override)
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

pfiles = dir('*-fio-fitparams.mat');
fioFit_total = [];

for j = 1:length(pfiles)
   infile = pfiles(j).name;
   load(infile,'fioFit');
   fioFit_total = [fioFit_total fioFit];
   clear fioFit;
end % (for j)

return;

