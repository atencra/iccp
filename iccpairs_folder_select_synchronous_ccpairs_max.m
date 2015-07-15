function iccpairs_folder_select_synchronous_ccpairs_max(basename)
%iccpairs_folder_select_synchronous_ccpairs_max Identify sharply synchronous
%neuron pairs in files
%
%   index = iccpairs_folder_select_synchronous_ccpairs(basename)
%
%   *-basename-isi-fixed-ccpairs.mat
%
%   basename indicates which ccpairs files to process. basename is a string
%   such as 'respcmb', 'resp', 'respcmbcmb', or 'strfcmb-pairs'
%
%   All file names that are processed will have isi-fixed-ccpairs.mat in
%   its name, since these files have been quality checked to ensure that
%   there are no inter-spike intervals < 0.75 ms.
% 
%   index = iccpairs_folder_select_synchronous_ccpairs(basename)


library('spike_train_cross_correlation');

s = sprintf('*-%s-isi-fixed-ccpairs.mat', basename);
d = dir(s);


for ii = 1:length(d)

    infile = d(ii).name;
    index = findstr(infile, '.mat');
    outfile = sprintf('%s-%s', infile(1:index-1), 'select-max.mat');
    fprintf('Infile: %s\n', infile);
    fprintf('Outfile: %s\n', outfile);

%     if ( ~exist(outfile,'file') )

        load(infile, 'ccpairs');

        index_max = [];
        index_chance = [];
        index_max_chance = [];

        for i = 1:length(ccpairs)

            exp = ccpairs(i).exp;
            site = ccpairs(i).site;
            chan = ccpairs(i).chan;
            model1 = ccpairs(i).model1;
            model2 = ccpairs(i).model2;
            position = ccpairs(i).position;
            fs = ccpairs(i).fs;

            pd = ccpairs(i).peakdelay;
            hw = ccpairs(i).halfwidth;
            ccc = ccpairs(i).ccc;

            conf_limit = ccpairs(i).conf_limit;
            delay = ccpairs(i).delay;
            q12 = ccpairs(i).q12;
            r12 = ccpairs(i).r12;

            spiketimes1 = ccpairs(i).spiketimes1;
            spiketimes2 = ccpairs(i).spiketimes2;

            nsp1 = length(spiketimes1);
            nsp2 = length(spiketimes2);


            % cross-correlation functions
            a = spiketimes1 / 1000;
            b = spiketimes2 / 1000;
            a = a(:)';
            b = b(:)';



            % Calculate if there are at least 100 spikes in each spike train
            if ( length(a) > 100 && length(b) > 100 )

                if ( length(a)>20000 ), a = a(1:20000); end
                if ( length(b)>20000 ), b = b(1:20000); end

                binwidth = 1.5; % bin resolution

                [ccdelay, ccraw, ccflat, ccnorm] = ...
                cc_crosscorr_binless(a, b, binwidth);

                index_ccdelay_end = find(ccdelay < -25 | ccdelay > 25);
                cc_chance = mean(ccnorm(index_ccdelay_end));

                % If the peak in the cross corr function is at 0, save the data
                [maxmax, imax] = max(ccnorm);
                imax = imax(1);
                if ( ccdelay(imax) == 0 )
                    index_max = [index_max i];
                    fprintf('Central peak.\n');
                end


                % Determine if central bin is 2 times greater than chance level
                index_zero = find(ccdelay == 0);
                cc_zero = ccnorm(index_zero);
                if ( cc_zero > 2 * cc_chance )
                    index_chance = [index_chance i];
                    fprintf('Central bin significant.\n');
                end


                % Determine if central bin is significant and if it is the maximum
                % in the cross corr function
                index_zero = find(ccdelay == 0);
                cc_zero = ccnorm(index_zero);
                if ( ( cc_zero > 2 * cc_chance ) && ccdelay(imax) == 0 )
                    index_max_chance = [index_max_chance i];
                    fprintf('Central peak and significant \n');
                end


                % Get spike waveforms
                % ----------------------------------------
                wavetime1 = ccpairs(i).meanwaveform1(:,1);
                waveform1 = ccpairs(i).meanwaveform1(:,2);
                waveform1= waveform1 ./ 2048;

                wavetime2 = ccpairs(i).meanwaveform2(:,1);
                waveform2 = ccpairs(i).meanwaveform2(:,2);
                waveform2= waveform2 ./ 2048;

                % Normalize the waveforms to [-1, 1]
                max1 = max([max(abs(waveform1)) max(abs(min(waveform1)))]);
                max2 = max([max(abs(waveform2)) max(abs(min(waveform2)))]);
                maxmax = max([max1 max2]);
                waveform1 = waveform1 ./ maxmax;
                waveform2 = waveform2 ./ maxmax;


                % Plot the waveforms
                % -------------------------------
                clf;
                subplot(2,1,1);
                hold on;
                plot(wavetime1, waveform1, 'k-', 'linewidth', 1.5);
                plot(wavetime2, waveform2, '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5);
                xlim([-1 2]);
                ylim([-1.05 1.05]);
                box off;
                tickpref;
                set(gca,'xtick', [-1:2], 'xticklabel', [-1:2]);
                xlabel('Time (ms)');
                set(gca,'yticklabel', [-1:1], 'ytick', [-1:1]);
                ylabel('Norm. Amp.');
                title(sprintf('#%.0f of %.0f: %s site%.0f chan%.0f model1[%.0f] model2[%.0f]', ...
                i, length(ccpairs), exp, site, chan, model1, model2));

                subplot(2,1,2);
                hold on;
                plot([0 0], [0 1.05*max(ccnorm)], '-', 'color', 0.6*ones(1,3), 'linewidth', 3);
                hb = bar(ccdelay, ccnorm);
                plot([-10 10], [cc_chance cc_chance], '-', 'color', 'r', 'linewidth', 3);
                plot([-10 10], 2*[cc_chance cc_chance], '--', 'color', 'r', 'linewidth', 3);
                set(hb,'facecolor', 0*ones(1,3));
                set(hb,'edgecolor', 0*ones(1,3));
                tickpref;
                xlim([-10 10]);
                xtick = -20:10:20;
                xtick = ccdelay;
                set(gca,'xtick', xtick, 'xticklabel', xtick);
                ylim([0 1.1*max(ccnorm)]);
                box off;
                xlabel('Delay (ms)');
                ylabel('CC');
                title(sprintf('Binless: Binwidth = %.2f ms', binwidth));

                set(gcf,'position', [150 141 400 700]);
                pause(0.1);

            end % (if)

        end % (for i)

        ccpairs_max = ccpairs(index_max);
        ccpairs_chance = ccpairs(index_chance);
        ccpairs_max_chance = ccpairs(index_max_chance);
        save(outfile, 'ccpairs_max', 'ccpairs_chance', 'ccpairs_max_chance');
        clear('ccpairs', 'ccpairs_max', 'ccpairs_chance', 'ccpairs_max_chance', ...
              'index_max', 'index_chance', 'index_max_chance');

%     else
%         fprintf('%s already processed.\n\n', outfile);
%     end

end % (for i = 1:length(d))


return;




