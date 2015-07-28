function iccp_batch_calc_fra_isi_fixed_crosscorr
% iccpairs_batch_calc_fra_isi_fixed_crosscorr Correlation between pairs of neurons for FRA stimulation
%
%     Processes all '*-fracmb-pairs-isi-fixed.mat' in the current
%     directory, and estimates the spike train cross-correlation between
%     all pairs of neurons on single electrode channels.
%
%     Data are saved in '*-fracmb-pairs-isi-fixed-ccpairs.mat' files
% 
%     iccpairs_batch_calc_fra_isi_fixed_crosscorr;

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

        [pairstrains] = iccp_get_spk_paired_data(spk);
        [ccpairs] = iccp_calc_spk_crosscorr(pairstrains);

        save(outfile, 'ccpairs');

        clear('s', 'spk', 'ccpairs');
        fprintf('Data saved in %s\n', outfile);

    else
        fprintf('%s already processed.\n\n', outfile);
    end

end % (for i)

return;


