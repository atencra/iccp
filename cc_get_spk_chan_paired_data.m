function [pairstrains] = cc_get_spk_chan_paired_data(spk)
% cc_get_spk_chan_paired_data(spk) Pairs from channels of spk
% 
%     [pairstrains] = cc_get_spk_chan_paired_data(spk)
%
%     spk is a struct array with the same number of elements. spk
%     holds the spike times.
% 
%     After spk is read in, for each channel of data, every
%     possible pair of neurons is computed, and the results are
%     saved in the struct array pairstrains. This struct holds the
%     spike trains for each pair of neurons.
% 
%
%     [pairstrains] = cc_get_spk_chan_paired_data(spk)


pairstrains = [];

cmb = nchoosek(1:length(spk),2); % determine all possible pairwise 
                                 % combinations for every recording
[nr,nc] = size(cmb);

for i = 1:nr

    index1 = cmb(i,1);
    index2 = cmb(i,2);

    % assign the data to a temporary struct variable
    datatemp.exp = spk(index1).exp;
    datatemp.site = spk(index1).site;
    datatemp.chan1 = spk(index1).chan;
    datatemp.chan2 = spk(index2).chan;

    datatemp.model1 = spk(index1).model;
    datatemp.model2 = spk(index2).model;

    datatemp.depth = spk(index1).depth;
    datatemp.position1 = spk(index1).position;
    datatemp.position2 = spk(index2).position;
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

end % (for i)


return;






