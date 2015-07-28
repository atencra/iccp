function [pairstrains] = iccp_get_spk_paired_data(spk)
% iccp_get_icc_fra_spk_paired_data Pairs from channels of spk
% 
%     [pairstrains] = iccp_get_icc_fra_spk_paired_data(spk);
% 
% 
%     Reads '*-tc1-fs*-spk-orig-isi-fixed.mat' files from the 
%     current directory. Each file holds spk, trigger, and fra data.
% 
%     spk is a struct array with the same number of elements. spk
%     holds the spike times.
% 
%     After spk is read in, for each channel of data, every
%     possible pair on a channel is computed, and the results are
%     saved in the struct array pairstrains.
% 
%
%     [pairstrains] = iccp_get_icc_fra_spk_paired_data(spk);


pairstrains = [];

chan = [spk.chan];
chan_unique = unique(chan);

for i = 1:length(chan_unique)

    index = find( chan_unique(i) == chan );

    cmb = nchoosek(index, 2); % determine all possible pairwise combinations

    [nr, nc] = size(cmb);

    if ( nr > 1 )

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

            if ( isfield(spk, 'meanwaveform') )
                datatemp.meanwaveform1 = spk(index1).meanwaveform;
                datatemp.meanwaveform2 = spk(index2).meanwaveform;
            else
                datatemp.meanwaveform1 = spk(index1).waveform;
                datatemp.meanwaveform2 = spk(index2).waveform;
            end

            datatemp.fs = spk(index1).fs;

            % assign the temporary struct to a larger struct array
            pairstrains = [pairstrains datatemp];

            clear('datatemp');

        end % (for j)
    end

end % (for i)


return;






