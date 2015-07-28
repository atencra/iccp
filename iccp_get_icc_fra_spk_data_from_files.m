function [pairstrains] = iccp_get_icc_fra_spk_data_from_files
% iccp_get_icc_fra_spk_data_from_files Pairs from channels of spk and fra data
% 
%     [pairstrains] = iccp_get_icc_fra_spk_data_from_files;
% 
%     Reads '*-tc1-fs*-spk-orig-isi-fixed.mat' files from the 
%     current directory. Each file holds spk, trigger, and fra data.
% 
%     spk and fra are struct array with the same number of elements. spk
%     holds the spike times, and fra holds the frequency response
%     areas.
% 
%     After spk and fra are read in, for each channel of data, every
%     possible pair on a channel is computed, and the results are
%     saved in the struct array pairstrains.
% 
%
%     [pairstrains] = iccp_get_icc_fra_spk_data_from_files





% For pairs of neurons from the same channel, get the data, and save it
% in a struct array for later processing
%
% The data that we'll save are the experiment, site, channel, model,
% depth, position, raw spike trains, and wave form durations, 

d = dir('*-tc1-fs*-spk-orig-isi-fixed.mat');

p = 0.01;

datastr = [];

for n = 1:length(d)

   filename = d(n).name;
   s = load(filename, 'spk');
   spk = s.spk; 
   chan = [spk.chan];
   chan_unique = unique(chan);


   index = findstr(filename, '-orig');
   prefix = filename(1:index-1);

   for i = 1:length(chan_unique)

      index = find( chan_unique(i) == chan );

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

         % assign the temporary struct to a larger struct array

         datastr = [datastr datatemp];
         clear('datatemp');

      end % (for j)

   end % (for i)

end % (for n)

return;






