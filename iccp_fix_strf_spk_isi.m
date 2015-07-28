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

        spknew = isi_filtered_spk(spk, delay);
        spk = spknew;

        save(outfile, 'strf', 'spk', 'trigger');
        clear('spknew', 'spk', 'strf', 'trigger');
        fprintf('Data saved in %s\n', outfile);

    else
        fprintf('%s already processed.\n\n', outfile);
    end

end %(for i)




function [center, range, nreps] = icc_exp_fra_params(exp, site, depth)
% icc_exp_fra_params TMS tuning curve parameters used in cat ICC experiments
% 
%     [center, range, nreps] = icc_exp_fra_params(exp, site, depth)
% 
%     Given the experiment date, site, and michigan probe insertion depth, 
%     this function returns the TMS system paramters that were used to 
%     present pure tones. 
% 
%     The TMS parameters are the center frequency of the tones presented,
%     the range of the frequencies, in octaves, and the the number
%     of tone repeats. For the TMS system, if a tone was repeated,
%     then the tone was played repeated in sequence. Thus, if tone A
%     was repeated five times, followed by tone B, then we would have
%     the stimulus sequence:
% 
%     A A A A A B B B B B ...
% 
%     center : center frequency, in kHz
% 
%     range : of frequencies, in octaves
% 
%     nreps : integer number of repeats, 1-N




if ( strcmp(exp, '2004-2-17') )

    if ( site == 1 && depth == 2395 )
        center = 1.6;
        range = 6;
        nreps = 5;

    elseif ( site == 4 && depth == 2405 )
        center = 1.6;
        range = 6;
        nreps = 5;

    elseif ( site == 11 && depth == 4046 )
        center = 2.5;
        range = 6;
        nreps = 5;


    elseif ( site == 12 && depth == 4270 )
        center = 5;
        range = 6;
        nreps = 5;

    elseif ( site == 13 && depth == 4057 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 14 && depth == 5000 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 14 && depth == 5000 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 19 && depth == 4767 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 19 && depth == 5005 )
        center = 2.5;
        range = 6;
        nreps = 3;
    else
        error('No FRA parameters for the exp/site/depth.');
    end

elseif ( strcmp(exp,'2003-11-24') )

    if ( ismember(site, 1:6) )
        center = 1.6;
        range = 6;
        nreps = 5;

    elseif ( ismember(site, 7:13) )
        center = 1;
        range = 6;
        nreps = 5;

    elseif ( site == 14 )
        if ( depth == 2420 )
            center = 1;
            range = 6;
            nreps = 5;
        elseif ( depth == 4620 )
        else
            error('Bad depth.');
        end
    else
        error('No FRA parameters for the exp/site/depth.');
    end

end



return;



