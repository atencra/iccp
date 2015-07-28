function iccp_fiofit_pairs(override)
% iccp_spktype_strf_crosscorr Correlation between pairs
%
% iccp_spktype_strf_crosscorr(override)
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


narginchk(0,1);
if ( nargin == 0 )
   override = 0;
end

pfiles = dir('*-strfcmb-pairs.mat');

for j = 1:length(pfiles)

   strffile = pfiles(j).name;

   fprintf('Getting STRF Params for %s\n', strffile);

   index = findstr(strffile, '.mat');
   basename = strffile(1:index-1);
   outfile = sprintf('%s-fiopairs.mat', basename );
   fprintf('Outfile = %s\n', outfile);

   d = dir(outfile);

   if ( length(d)==0 || override ) % We haven't previously processed the data ...

      load(strffile, 'spk');

      rtfparamsfile = sprintf('%s-sta-fio.mat', basename );
      load(rtfparamsfile, 'fio');

      rtfparamsfile = sprintf('%s-sta-fio-fitparams.mat', basename );
      load(rtfparamsfile, 'fioFit');


      % check to make sure that each struct has the same number of elements
      nFiles = length(spk);
      totalFiles = length(spk) + length(fio) + length(fioFit);

      if ( 3*nFiles ~= totalFiles )
         error('Mismatch between number of files for different analyses.');
      end

      [fiopairs] = get_fiopairs_spiketrains(spk, fio, fioFit);

      fprintf('Saving Fio Pairs in %s\n\n', outfile);
      save(outfile, 'fiopairs');
      clear('spk', 'fio', 'fioFit');

   else
       fprintf('Outfile = %s already exists.\n', outfile);
   end % (if length(d)==0)

end % (for j)

return;



%------------------- Function Definitions ------------------------

function [datastr] = get_fiopairs_spiketrains(spk, fio, fioFit)
% For pairs of neurons from the same channel, get the data, and save it
% in a struct array for later processing
%
% The data that we'll save are the experiment, site, channel, model,
% depth, position, raw spike trains, and wave form durations,


p = 0.01;

datastr = [];

chan = [spk.chan];
chan_unique = unique(chan);

for i = 1:length(chan_unique)

   index = find( chan_unique(i) == chan );

   if ( index > 1 )

      cmb = nchoosek(index, 2); % determine all possible pairwise combinations

      [nr, nc] = size(cmb);

      for j = 1:nr

         index1 = cmb(j,1);
         index2 = cmb(j,2);


         % assign the data to a temporary struct variable
         datatemp.exp = spk(index1).exp;
         datatemp.site = spk(index1).site;
         datatemp.chan = spk(index1).chan;

         datatemp.model1 = spk(index1).model;
         datatemp.model2 = spk(index2).model;

         datatemp.depth = spk(index1).depth;
         datatemp.position = spk(index1).position;
         datatemp.stim = spk(index1).stim;
         datatemp.atten = spk(index1).atten;

         datatemp.header1 = spk(index1).header;
         datatemp.header2 = spk(index2).header;

         datatemp.waveform1 = spk(index1).waveform;
         datatemp.waveform2 = spk(index2).waveform;

         datatemp.spiketimes1 = spk(index1).spiketimes;
         datatemp.spiketimes2 = spk(index2).spiketimes;

         datatemp.outliers1 = spk(index1).outliers;
         datatemp.outliers2 = spk(index2).outliers;

         datatemp.meanwaveform1 = spk(index1).meanwaveform;
         datatemp.meanwaveform2 = spk(index2).meanwaveform;

         datatemp.fs = spk(index1).fs;


         % Get the nonlinearity parameters
         fio1 = fio(index1);
         fio2 = fio(index2);

         [sta1Mn, pspkx1Mn, x1Mn, y1max] = get_pairs_sta_fio_stats(fio1);
         [sta2Mn, pspkx2Mn, x2Mn, y2max] = get_pairs_sta_fio_stats(fio2);

         [x12, i1, i2] = intersect(x1Mn, x2Mn);

         pspkx1 = pspkx1Mn(i1);
         pspkx2 = pspkx2Mn(i2);

         r = corrcoef(pspkx1, pspkx2); % correlation between nonlinearities
         fiosimilarity = r(1,2);

         asi1 = fio_asi(x12, pspkx1);
         asi2 = fio_asi(x12, pspkx2);

         datatemp.fiosimilarity = fiosimilarity;
         datatemp.fioasi1 = asi1;
         datatemp.fioasi2 = asi2;


         fiofit1 = fioFit(index1);
         fiofit2 = fioFit(index2);

         theta1 = fiofit1.fitParams_sta(2);
         sigma1 = fiofit1.fitParams_sta(3);
         nmse1 = fiofit1.nmse_sta;
         r21 = fiofit1.r2_sta;

         theta2 = fiofit2.fitParams_sta(2);
         sigma2 = fiofit2.fitParams_sta(3);
         nmse2 = fiofit2.nmse_sta;
         r22 = fiofit2.r2_sta;

         datatemp.theta1 = theta1;
         datatemp.sigma1 = sigma1;
         datatemp.nmse1 = nmse1;
         datatemp.r21 = r21;

         datatemp.theta2 = theta2;
         datatemp.sigma2 = sigma2;
         datatemp.nmse2 = nmse2;
         datatemp.r22 = r22;


         % assign the temporary struct to a larger struct array

         datastr = [datastr datatemp];
         clear('datatemp');

      end % (for j)

   end % (if)

end % (for i)

return;



