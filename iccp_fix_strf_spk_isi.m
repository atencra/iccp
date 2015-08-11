function iccp_fix_strf_spk_isi
% iccp_fix_strf_spk_isi Modify spike train files so no ISIs < 0.5 ms
% 
%     iccp_fix_strf_spk_isi
% 
%     Fixes spike train data saved in spk struct array data. Spike trains
%     are examined, and if two spikes are separated by less than 0.5 ms, 
%     the second spike is removed from the spike train.
% 
%     The function looks for *-strfcmb-pairs.mat files in the current
%     directory and fixes the data inside them. 
% 
%     New data files will be saved in the form 
% 
%         *-strfcmb-pairs-isi-fixed.mat
% 
%     where the ending -isi-fixed.mat has been appended to the 
%     old file name.
% 
%     iccp_fix_strf_spk_isi;


library('spikesort');

d = dir( sprintf('*-strfcmb-pairs.mat') );

delay = 0.75;

for i = 1:length(d)

    fprintf('\n#%.0f of %.0f\n', i, length(d));
    fprintf('%s\n', d(i).name);

    file = d(i).name;
    index = findstr(file, '.mat');
    outfile = [file(1:index-1) '-isi-fixed.mat'];
    fprintf('%s\n', outfile);

    if ( ~exist(outfile, 'file') )

        s = load(file, 'strf','spk', 'trigger');
        spk = s.spk;
        strf = s.strf;
        trigger = s.trigger;

        spknew = iccp_isi_filtered_spk(spk, delay);
        spk = spknew;

        save(outfile, 'strf', 'spk', 'trigger');
        clear('spknew', 'spk', 'strf', 'trigger');
        fprintf('Data saved in %s\n', outfile);

    else
        fprintf('%s already processed.\n\n', outfile);
    end

end %(for i)

return;







