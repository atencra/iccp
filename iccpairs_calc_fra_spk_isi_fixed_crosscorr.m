function iccpairs_calc_fra_spk_isi_fixed_crosscorr
% iccpairs_calc_fra_spk_isi_fixed_crosscorr Correlation between pairs/cell types
%
% 
% iccpairs_calc_fra_spk_isi_fixed_crosscorr
% --------------------------------------------------------
%
% The function finds the neurons that were recorded from the 
% same channel, and then calculates the cross correlations
% between the neurons.
%
% The function returns three celltype based struct arrays: 
% fsufsu - holds the fsu-fsu pair spike trains; fsursu - holds the 
% fsu-rsu pair spike trains; and rsursu hold the rsu-rsu pair data
%
% The fourth returned argument, datastr, holds all the pair spike trains, 
% before the data has been separated in celltype.
%
% You may also input these variables after you the function the first time
% to save compuation time.
%
% The function looks for files in the current directory that have the 
% following form:
% 
%        *-pairs-orig-isi-fixed.mat
%
% caa 2/4/12
%
%     iccpairs_calc_fra_spk_isi_fixed_crosscorr


narginchk(0,0);

d = dir('*-tc1-fs*-spk-orig-isi-fixed.mat');

for i = 5:length(d)

    fprintf('File %.0f of %.0f\n', i, length(d));

    filename = d(i).name;
    index = findstr(filename, '.mat');
    basename = filename(1:index-1);
    outfile = sprintf('%s-ccpairs.mat',basename)

    s = load(filename, 'spk');
    spk = s.spk; 
    chan = [spk.chan];
    chan_unique = unique(chan);

    [pairstrains] = iccpairs_get_icc_fra_spk_paired_data(spk);

    [ccpairs] = iccpairs_calc_icc_fra_spk_crosscorr(pairstrains);

%     save(outfile, 'ccpairs');

end % (for i)


return;





