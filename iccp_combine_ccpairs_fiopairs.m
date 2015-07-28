function newpairs = iccp_combine_ccpairs_fiopairs(ccpairs, fiopairs)
% iccp_combine_ccpairs_fiopairs Correlation between pairs
%
% iccp_combine_ccpairs_fiopairs(ccpairs, fiopairs)
% -----------------------------------------------------------
%
% Reads in data from pairs analysis and estimates the cross-correlation
% between the spike trains of each neuron in a pair. Data is read from 
% files in a directory that have already been processed. The data files
% are:
%
% *-strfcmb-pairs.mat, 
% *-strfcmb-pairs-params.mat,
% *-strfcmb-pairs-ptparams.mat, 
% *-strfcmb-pairs-rtf-params.mat,
% *-strfcmb-pairs-sta-fio.mat
%
% Each file contains a struct array. The number of struct elements in each
% struct array for each file must match.
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


narginchk(2,2);

if ( length(ccpairs) ~= length(fiopairs) )
   error('ccpairs and fiopairs must be the same length.');
end

newpairs = ccpairs;


for i = 1:length(ccpairs)

   if ( strcmp(ccpairs(i).exp, fiopairs(i).exp) && ccpairs(i).site == fiopairs(i).site && ...
      ccpairs(i).chan == fiopairs(i).chan && sum(ccpairs(i).model1 - fiopairs(i).model1) == 0 && ...
      sum(ccpairs(i).model2 - fiopairs(i).model2) == 0 && strcmp(ccpairs(i).stim, fiopairs(i).stim) )
      
      newpairs(i).theta1 = fiopairs(i).theta1;
      newpairs(i).sigma1 = fiopairs(i).sigma1;
      newpairs(i).nmse1 = fiopairs(i).nmse1;
      newpairs(i).r21 = fiopairs(i).r21;

      newpairs(i).theta2 = fiopairs(i).theta2;
      newpairs(i).sigma2 = fiopairs(i).sigma2;
      newpairs(i).nmse2 = fiopairs(i).nmse2;
      newpairs(i).r22 = fiopairs(i).r22;

   else
      error('ccpairs and fiopairs do not match.');
   end

end % (for i)


















