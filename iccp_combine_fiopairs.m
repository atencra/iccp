function fiopairs_total = iccp_combine_fiopairs
% iccp_combine_fiopairs Gather FIO nonlinearity data into one struct array
%
% fiopairs_total = iccp_combine_fiopairs;
% -----------------------------------------------------------
%
% Simple function to gather all the *-strfcmb-pairs-fiopair.mat data in the
% current directory and append it to a larger struct array database.


narginchk(0,0);

pfiles = dir('*-strfcmb-pairs-fiopairs.mat');
fiopairs_total = [];

for j = 1:length(pfiles)
   fiofile = pfiles(j).name;
   load(fiofile,'fiopairs');
   fiopairs_total = [fiopairs_total fiopairs];
   clear fiopairs;
end % (for j)

return;

