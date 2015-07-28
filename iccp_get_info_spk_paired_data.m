function [pairstrains] = iccp_get_info_spk_paired_data(resp)
% iccp_get_icc_fra_spk_paired_data Pairs from channels of spk
% 
%     [pairstrains] = iccp_get_info_spk_paired_data(resp)
% 
%     Returns pairs of spike trains that are stored in the struct array
%     resp. Every pairwise combination of neurons in resp is calculate.
%     The combination are calculated for each pair of neurons on the 
%     same recording channel.
%
%     resp holds the spike times in response to a repeated dmr stimulus 
%     sequence.
%
%
%     [pairstrains] = iccp_get_info_spk_paired_data(resp);


pairstrains = [];

chan = [resp.chan];
chan_unique = unique(chan);

for i = 1:length(chan_unique)

    index = find( chan_unique(i) == chan );

    if ( length(index) > 1 )
        cmb = nchoosek(index, 2); % determine all possible pairwise combinations
        [nr, nc] = size(cmb);
    else
        nr = 0;
    end

    if ( nr > 1 )

        for j = 1:nr

            index1 = cmb(j,1);
            index2 = cmb(j,2);

            % assign the data to a temporary struct variable
            datatemp.exp = resp(index1).exp;
            datatemp.site = resp(index1).site;
            datatemp.chan = resp(index1).chan;

            datatemp.model1 = resp(index1).model;
            datatemp.model2 = resp(index2).model;

            datatemp.depth = resp(index1).depth;
            datatemp.position = resp(index1).position;
            datatemp.stim = resp(index1).stim;
            datatemp.atten = resp(index1).atten;

            datatemp.header1 = resp(index1).header;
            datatemp.header2 = resp(index2).header;

            datatemp.waveform1 = resp(index1).waveform;
            datatemp.waveform2 = resp(index2).waveform;

            datatemp.spiketimes1 = resp(index1).spiketimes;
            datatemp.spiketimes2 = resp(index2).spiketimes;

            datatemp.outliers1 = resp(index1).outliers;
            datatemp.outliers2 = resp(index2).outliers;

            if ( isfield(resp, 'meanwaveform') )
                datatemp.meanwaveform1 = resp(index1).meanwaveform;
                datatemp.meanwaveform2 = resp(index2).meanwaveform;
            else
                datatemp.meanwaveform1 = resp(index1).waveform;
                datatemp.meanwaveform2 = resp(index2).waveform;
            end

            datatemp.fs = resp(index1).fs;

            % assign the temporary struct to a larger struct array
            pairstrains = [pairstrains datatemp];

            clear('datatemp');

        end % (for j)
    end

end % (for i)


return;






