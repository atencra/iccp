function iccp_fix_info_spk_isi
% iccp_fix_info_spk_isi Modify spike train files so no ISIs < 0.75 ms
% 
%     iccpairs_fix_info_spk_isi
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

datatype = {'resp', 'respcmb', 'respcmbcmb'};
delay = 0.75;

for j = 1:length(datatype)

    d = dir( sprintf('*-dmrrep-fs18115-%s.mat',datatype{j}) );

    for i = 1:length(d)

        fprintf('\n#%.0f of %.0f\n', i, length(d));
        fprintf('%s\n', d(i).name);

        file = d(i).name;
        index = findstr(file, '.mat');
        outfile = [file(1:index-1) '-isi-fixed.mat'];
        fprintf('%s\n', outfile);

        if ( ~exist(outfile, 'file') )

            s = load(file, 'resp');
            resp = s.resp;

            for ii = 1:length(resp)
                spiketimes = resp(ii).spiketimes;
                newspiketimes = isi_filtered_spiketimes(spiketimes, delay);
                resp(ii).spiketimes = newspiketimes;
            end % (for ii)

            save(outfile, 'resp');
            clear('resp');
            fprintf('Data saved in %s\n', outfile);

        else
            fprintf('%s already processed.\n\n', outfile);
        end

    end %(for i)

end %(for j)


return;







