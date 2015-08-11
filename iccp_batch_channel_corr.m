function iccp_batch_channel_corr
% iccp_batch_channel_corr Batch compute correlation between pairs of neurons
%
%     iccp_batch_channel_corr
%     --------------------------------------------------------
% 
%     The function finds the neurons that were recorded from the 
%     same channel, and then calculates the cross correlations
%     between the neurons.
% 
%     The data that are used are contained in:
% 
%             '*-strfcmb-pairs-isi-fixed.mat'
% 
%     The spk struct array data in these files has been checked to 
%     make sure that no spike train has a inter-spike interval less than 0.5 ms.
% 
%     This makes the cross-correlation results more interpretable,
%     since peaks near 0 can only be a result of synchronous spiking,
%     and not spike-sorting errors.
% 
%     The function processes all '*-strfcmb-pairs-isi-fixed.mat' files
%     using the function calls:
% 
%             [pairstrains] = iccp_get_spk_paired_data(spk);
%             [ccpairs] = iccp_calc_spk_crosscorr(pairstrains);
% 
%     The data are saved in files of the form:
% 
%             '*-strfcmb-pairs-isi-fixed-ccpairs.mat'
% 
% 
%     iccpairs_calc_isi_fixed_crosscorr;

narginchk(0,0);

d = dir('*-strfcmb-pairs-isi-fixed.mat');

for i = 1:length(d)

    fprintf('\nFile %.0f of %.0f\n', i, length(d));

    filename = d(i).name;
    index = findstr(filename, '.mat');
    basename = filename(1:index-1);
    outfile = sprintf('%s-ccpairs.mat',basename);
    fprintf('%s\n', outfile);

   if ( ~exist(outfile, 'file') )

        str = load(filename, 'spk');

        if ( length(str.spk) > 1 )

            [pairstrains] = iccp_get_paired_spk_data(str.spk);
            [ccpairs] = iccp_calc_spk_crosscorr(pairstrains);

            save(outfile, 'ccpairs');

            clear('s', 'spk', 'pairstrains', 'ccpairs');
            fprintf('Data saved in %s\n', outfile);

        end 

    else
        fprintf('%s already processed.\n\n', outfile);
    end

end % (for i)

return;


