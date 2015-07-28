function iccp_folder_select_synchronous_ccpairs(basename)
%iccp_folder_select_synchronous_ccpairs Identify sharply synchronous
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
%   iccp_folder_select_synchronous_ccpairs(basename)


library('agmon');

s = sprintf('*-%s-isi-fixed-ccpairs.mat', basename);
d = dir(s);

for ii = 1:length(d)

    infile = d(ii).name;
    load(infile, 'ccpairs');
    index = findstr(infile, '.mat');
    outfile = sprintf('%s-%s', infile(1:index-1), 'select.mat');
    fprintf('Infile: %s\n', infile);
    fprintf('Outfile: %s\n', outfile);

    if ( ~exist(outfile,'file') )

        figure;

        index = [];
        index_max = [];

        for i = 1:length(ccpairs)

            clf;

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
            if ( length(a)>20000 ), a = a(1:20000); end
            if ( length(b)>20000 ), b = b(1:20000); end
            binwidth = 1.5; % bin resolution
            [ccdelay, ccraw, ccflat, ccnorm] = ...
                agmon_crosscorr_binless(a, b, binwidth);

            % autocorrelation functions
            acbinwidth = 0.5; % bin resolution
            [accdelay, accraw, accflat, accnorm] = ...
                agmon_crosscorr_binless(a, a, acbinwidth);

            [bccdelay, bccraw, bccflat, bccnorm] = ...
                agmon_crosscorr_binless(b, b, acbinwidth);


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
           subplot(3,2,1);
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


           % Plot the cross-covariance
           % -------------------------------
            subplot(3,2,2);
            upper95qab = conf_limit;
            lower95qab = -conf_limit;
            hold on;
            ymin = min([min(q12) lower95qab]);
            ymax = max([max(q12) upper95qab]);
            range = ymax-ymin;
            plot([0 0], [ymin-0.1*range ymax+0.1*range], '-', ...
                'linewidth', 3, 'color', 0.5*ones(1,3));
            plot([min(delay) max(delay)], [0 0], '-', ...
                'linewidth', 3, 'color', 0.5*ones(1,3));
            bar(delay, q12, 'k');

            plot([min(delay) max(delay)], [upper95qab upper95qab], 'r-');
            plot([min(delay) max(delay)], [lower95qab lower95qab], 'r-');

            plot([pd pd], [ymin-0.1*range ymax+0.1*range], 'r-');

            xlim([min(delay) max(delay)]);
            xlim([-25 25]);
            xtick = -20:10:20;
            set(gca,'xtick', xtick, 'xticklabel', xtick);
            xlabel('Delay (ms)');
            ylabel('Cross Covariance');
            ylim([ymin-0.1*range ymax+0.1*range]);
            tickpref;
            title(sprintf('#%.0f of %.0f: PD = %.2f, HW = %.2f', ...
                i, length(ccpairs), pd, hw));



            subplot(3,2,4);
            hold on;
            plot([0 0], [0 1.1*max(r12)], '-', ...
                'linewidth', 3, 'color', 0.5*ones(1,3));
            bar(delay, r12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
            plot([pd pd], [0 1.05*max(r12)], 'r-');
            xlim([min(delay) max(delay)]);
            xlim([-25 25]);
            ylim([0 1.1*max(r12)]);
            xtick = -20:10:20;
            set(gca,'xtick', xtick, 'xticklabel', xtick);
            xlabel('Delay (ms)');
            ylabel('Cross-correlation');
            %    ylim([ymin-0.1*range ymax+0.1*range]);
            tickpref;


            subplot(3,2,6);
            hold on;
            plot([0 0], [0 1.05*max(ccnorm)], '-', 'color', 0.6*ones(1,3), 'linewidth', 3);
            hb = bar(ccdelay, ccnorm);
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

            [maxmax, imax] = max(ccnorm);
            imax = imax(1);
            if ( ccdelay(imax) == 0 )
                index_max = [index_max i];
            end



            subplot(3,2,3);
            hold on;
            plot([0 0], [0 1.05*max(accnorm)], '-', 'color', 0.6*ones(1,3), 'linewidth', 3);
            hb = bar(accdelay, accnorm);
            set(hb,'facecolor', 0*ones(1,3));
            set(hb,'edgecolor', 0*ones(1,3));
            tickpref;
            xlim([-3 3]);
            xtick = accdelay;
            set(gca,'xtick', xtick, 'xticklabel', xtick);
            ylim([0 1.1*max(accnorm)]);
            box off;
            xlabel('Delay (ms)');
            ylabel('Auto Corr A');

            subplot(3,2,5);
            hold on;
            plot([0 0], [0 1.05*max(bccnorm)], '-', 'color', 0.6*ones(1,3), 'linewidth', 3);
            hb = bar(bccdelay, bccnorm);
            set(hb,'facecolor', 0*ones(1,3));
            set(hb,'edgecolor', 0*ones(1,3));
            tickpref;
            xlim([-3 3]);
            xtick = bccdelay;
            set(gca,'xtick', xtick, 'xticklabel', xtick);
            ylim([0 1.1*max(bccnorm)]);
            box off;
            xlabel('Delay (ms)');
            ylabel('Auto Corr B');

            set(gcf,'position', [150 141 700 700]);

            fprintf('\n');
            inp = input(sprintf('#%.0f of %.0f: Keep or reject? <enter> = reject, 1 = keep: ',...
                i,length(ccpairs)));
            if ( ~isempty(inp) )
                index = [index i];
                save('temp-ccpairs-select', 'index');
            end

            clf;

        end % (for i)

        ccpairs_select = ccpairs(index); % rename and save
        ccpairs_max = ccpairs(index_max); 
        save(outfile, 'ccpairs_select', 'ccpairs_max');
        clear('ccpairs');
        fprintf('Data saved in: %s\n\n', outfile);
    else
        fprintf('%s already processed.\n\n', outfile);
    end

end % (for i = 1:length(d))


return;




