function iccpairs_batch_calc_fra_isi_fixed_crosscorr
% iccpairs_batch_calc_fra_isi_fixed_crosscorr Correlation between pairs of neurons
%
%     iccpairs_batch_calc_fra_isi_fixed_crosscorr
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
%     The spk struct array data in these files has been processed so
%     that no spike train has a inter-spike interval less than 0.5 ms.
% 
%     This will make the cross-correlation results more interpretable,
%     since peaks near 0 can only be a result of synchronous spiking,
%     and not spike-sorting errors.
% 
%     The function processes all '*-strfcmb-pairs-isi-fixed.mat' files
%     using the function calls:
% 
%             [pairstrains] = iccpairs_get_spk_paired_data(spk);
%             [ccpairs] = iccpairs_calc_icc_fra_spk_crosscorr(pairstrains);
% 
%     The data are saved in files of the form:
% 
%             '*-strfcmb-pairs-isi-fixed-ccpairs.mat'
% 
% 
%     iccpairs_calc_isi_fixed_crosscorr;

narginchk(0,0);

d = dir('*-fracmb-pairs-isi-fixed.mat');

for i = 1:length(d)

    fprintf('\nFile %.0f of %.0f\n', i, length(d));

    infile = d(i).name;
    index = findstr(infile, '.mat');
    basename = infile(1:index-1);
    outfile = sprintf('%s-ccpairs.mat',basename);
    fprintf('Infile:  %s\n', infile);
    fprintf('Outfile: %s\n', outfile);

    if ( ~exist(outfile, 'file') )

        s = load(infile, 'spk');
        spk = s.spk; 

        [pairstrains] = iccpairs_get_spk_paired_data(spk);
        [ccpairs] = iccpairs_calc_spk_crosscorr(pairstrains);

        save(outfile, 'ccpairs');

        clear('s', 'spk', 'ccpairs');
        fprintf('Data saved in %s\n', outfile);

    else
        fprintf('%s already processed.\n\n', outfile);
    end

end % (for i)

return;


