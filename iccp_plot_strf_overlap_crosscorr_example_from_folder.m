function iccp_plot_strf_overlap_crosscorr_example_from_folder(ccoverlap)
% iccp_plot_strf_overlap_crosscorr_example_from_folder STRF overlap for two neurons
% 
%     iccp_plot_strf_overlap_crosscorr_example(ccoverlap)
%     ------------------------------------------------------------------
%     Plot example comparison between cross-covariance functions and the
%     correlation between STRF slices.
%
%     The current folder is searched for the STRFs that match the current
% 
%     strf : struct array of strfs, where some channels have two STRFs.
% 
%     trigger : vector of trigger times, in sample number.
% 
%     ccpairs : struct array holding pairwise correlation data for the
%     channels in STRF. ccpairs holds the cross-covariance function, the
%     correlogram, the peak delay, half-width, etc. from correlation
%     analysis.
% 
% 
%     The function finds channels where there are at least two STRFs, and
%     then correlates the STRFs using a slice through the CF. Each pairwise
%     combination of neurons from the same channel in strf is plotted.
% 
%     The RF slice correlation is then compared to the correlation half-
%     width to see if the STRF overlap predicts the spike train correlation
%     width.



% faxis = strf(1).faxis;
% taxis = strf(1).taxis;
% dtemp = taxis(2) - taxis(1);
% if ( dtemp < 0.1 )
%     taxis = taxis * 1000; % convert to ms
% end
% 
% dt = ( taxis(2) - taxis(1) ); % size of bins
% maxlag = ceil(50 / dt); % Number of bins to cover 50 ms



% Earlier code:
% % Part 1. Get slice through peaks:
% 
% [i1, j1] =  find(max(max(strf1)) == strf1);
% slice1 = strf1(i1,:);
% slice11 = slice1 ./ max([abs(slice1) eps]); % strf1, peak1
% 
% [i2, j2] =  find(max(max(strf2)) == strf2);
% slice2 = strf2(i2,:);
% slice22 = slice2 ./ max([abs(slice2) eps]); % strf2, peak 2
% 
% 
% % slice2 = strf2(i1,:);
% % slice21 = slice2 ./ max([abs(slice2) eps]); % strf 2, peak 1
% % 
% % slice1 = strf1(i2,:)
% % slice12 = slice1 ./ max([abs(slice1) eps]); % strf 1, peak 2
% 
% % [slice11(:) slice22(:) slice21(:) slice12(:)]
% 
% 
% % [c12, lags] = xcorr(slice11, slice22, maxlag, 'coeff');
% % [c11, lags] = xcorr(slice11, slice21, maxlag, 'coeff');
% % [c22, lags] = xcorr(slice12, slice22, maxlag, 'coeff');
% 
% [c12, lags] = xcorr(slice11, slice22, maxlag);
% % [c11, lags] = xcorr(slice11, slice21, maxlag);
% % [c22, lags] = xcorr(slice12, slice22, maxlag);
% c12 = c12 ./ max(c12);
% % c11 = c11 ./ max(c11);
% % c22 = c22 ./ max(c22);
% 
% 
% lags = lags * dt;
% 
% 
% index12 = find(c12 == max(c12));
% peak12 = c12(index12);
% lag_peak12 = lags(index12);
% 
% index_low12 = find(c12 < 0.5 *peak12 & lags < lag_peak12 );
% index_low12 = max([1 index_low12]);
% lag_low12 = lags( max( index_low12 ) );
% 
% index_high12 = find(c12 < 0.5 *peak12 & lags > lag_peak12 );
% index_high12 = min([length(c12) index_high12]);
% lag_high12 = lags( min( index_high12 ) );
% 
% 
% 
% 
% clf;
% 
% % Part 1.
% subplot(4,3,1);
% plot_strf_symmetric_colormap(strf1);
% title('STRF1');
% 
% subplot(4,3,2);
% plot_strf_symmetric_colormap(strf2);
% title('STRF2');
% 
% subplot(4,3,4);
% plot(taxis, slice11, 'k-');
% xlim([min(taxis) max(taxis)]);
% ylim([min(slice11) max(slice11)]);
% box on;
% tickpref;
% title('STRF1 peak, STRF1 slice');
% 
% subplot(4,3,5);
% plot(taxis, slice22, 'k-');
% xlim([min(taxis) max(taxis)]);
% ylim([min(slice22) max(slice22)]);
% box on;
% tickpref;
% title('STRF2 peak, STRF2 slice');
% 
% subplot(4,3,6);
% hold on;
% patch([lag_low12 lag_low12 lag_high12 lag_high12],1*[-1 1 1 -1],0.85*[1 1 1]);
% patch(1.25*[-1 -1 1 1],1*[-1 1 1 -1],[0.6 0.6 0.6]);
% plot(lags,c12,'k-');
% xtick = sort([-50 -25 -10 10 25 50]);
% set(gca,'xtick', xtick, 'xticklabel', xtick);
% xlim([min(lags) max(lags)]);
% box on;
% tickpref;
% title('Slice through each peak');


for ii = 1:length(ccoverlap)

    exp = ccoverlap(ii).exp;
    site = ccoverlap(ii).site;
    chan = ccoverlap(ii).chan;
    depth = ccoverlap(ii).depth;
    atten = ccoverlap(ii).atten;
    stim = ccoverlap(ii).stim;
    model1 = ccoverlap(ii).model1;
    model2 = ccoverlap(ii).model2;

    filePattern = sprintf('%s-site%.0f-%.0fum-%.0fdb-%s-fs*-strfcmb-pairs-isi-fixed.mat',...
                    exp, site, depth, atten, stim);

    strfFile = dir(filePattern);

    if ( length(strfFile) == 1 )

        load(strfFile.name);
        fprintf('%s\n', strfFile.name);


        % Get the STRFs for the ccoverlap element
        [strf1, strf2] = pairs_get_strf_pair_for_crosscorr(strf, ccoverlap(ii));

        if ( ~isempty(strf1) && ~isempty(strf2) )

            % Now get the STRF similarity data
            taxis = strf1.taxis;
            faxis = strf1.faxis;

            rf1 = strf1.rfcontra;
            rf2 = strf2.rfcontra;

            n01 = strf1.n0contra;
            n02 = strf2.n0contra;

            mdb = strf1.mdb;
            fs = strf1.fs;
            stimdur = ( trigger(end) - trigger(1) ) / fs;
            soundtype = strf1.stim;
            p = 0.01;

            rfsig1 = significant_strf(rf1, p, n01, mdb, stimdur, soundtype);
            rfsig2 = significant_strf(rf2, p, n02, mdb, stimdur, soundtype);

            [overlap] = strf_overlap(rfsig1, rfsig2, taxis);

            dtemp = taxis(2) - taxis(1);
            if ( dtemp < 0.1 )
                taxis = taxis * 1000; % convert to ms
            end
            dt = ( taxis(2) - taxis(1) ); % size of bins
            maxlag = ceil(25 / dt); % Number of bins to cover 50 ms


            
            indt = find(taxis >= -0 & taxis <= 25);
            taxispos = taxis(indt);
            rfsig1 = rfsig1(:,indt);
            rfsig2 = rfsig2(:,indt);

            indf = find(faxis >= 500 & faxis <= 20000);
            faxispos = faxis(indf) / 1000;
            rfsig1 = rfsig1(indf,:);
            rfsig2 = rfsig2(indf,:);

            cc_delay = ccoverlap(ii).delay;
            cc_q12 = ccoverlap(ii).q12;
            cc_q12 = cc_q12 ./ max(cc_q12);


            % Correlate the receptive fields
            [i1, j1] =  find(max(max(rfsig1)) == rfsig1);
            slice1 = rfsig1(i1,:);
            slice11 = slice1 ./ max([abs(slice1) eps]); % strf1, peak1

            [i2, j2] =  find(max(max(rfsig2)) == rfsig2);
            slice2 = rfsig2(i2,:);
            slice22 = slice2 ./ max([abs(slice2) eps]); % strf2, peak 2


            [c12, lags] = xcorr(slice11, slice22, maxlag);
            c12 = c12 ./ max(c12);

            lags = lags * dt;

            index12 = find(c12 == max(c12));
            peak12 = c12(index12);
            lag_peak12 = lags(index12);

            index_low12 = find(c12 < 0.5 *peak12 & lags < lag_peak12 );
            index_low12 = max([1 index_low12]);
            lag_low12 = lags( max( index_low12 ) );

            index_high12 = find(c12 < 0.5 *peak12 & lags > lag_peak12 );
            index_high12 = min([length(c12) index_high12]);
            lag_high12 = lags( min( index_high12 ) );





            cmap = brewmaps('rdbu', 15);


            figure;

            % Part 1.
            subplot(5,1,1);
            hold on;
            minmin = min(min(rfsig1));
            maxmax = max(max(rfsig1));
            boundary = 0.9 * max([abs(minmin) abs(maxmax)]);
            imagesc(rfsig1);
            xtick = [1 ceil(length(taxispos)/2) length(taxispos)];
            xticklabel = taxispos(xtick);
            xticklabel = ceil(xticklabel*10)/10;
            set(gca,'xtick', xtick, 'xticklabel', xticklabel);
            xlim([1 length(taxispos)]);

            ytick = [1 ceil(length(faxispos)/2) length(faxispos)];
            yticklabel = faxispos(ytick);
            yticklabel = round(yticklabel*10)/10;
            set(gca,'ytick', ytick, 'yticklabel', yticklabel);
            ylim([1 length(faxispos)]);


            ylim([1 length(faxispos)]); 
            set(gca,'ydir', 'normal');
            set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
            set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
            colormap(cmap);
            plot([xlim], [i1 i1], 'k-');


            subplot(5,1,2);
            hold on;
            minmin = min(min(rfsig2));
            maxmax = max(max(rfsig2));
            boundary = 0.9 * max([abs(minmin) abs(maxmax)]);
            imagesc(rfsig2);

            xtick = [1 ceil(length(taxispos)/2) length(taxispos)];
            xticklabel = taxispos(xtick);
            xticklabel = ceil(xticklabel*10)/10;
            set(gca,'xtick', xtick, 'xticklabel', xticklabel);
            xlim([1 length(taxispos)]);

            ytick = [1 ceil(length(faxispos)/2) length(faxispos)];
            yticklabel = faxispos(ytick);
            yticklabel = round(yticklabel*10)/10;
            set(gca,'ytick', ytick, 'yticklabel', yticklabel);
            ylim([1 length(faxispos)]);


            ylim([1 length(faxispos)]); 
            set(gca,'ydir', 'normal');
            set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
            set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
            colormap(cmap);
            plot([xlim], [i2 i2], 'k-');




            subplot(5,1,3);
            plot(taxispos, slice11, 'k-');
            xlim([min(taxispos) max(taxispos)]);
            ylim([min(slice11) max(slice11)]);
            xtick = [1 ceil(length(taxispos)/2) length(taxispos)];
            xtick = taxispos(xtick);
            xticklabel = ceil(xtick*10)/10;
            set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
            box on;
            tickpref;
            ylabel('Norm. Amplitude');
            title('STRF1 slice');


            subplot(5,1,4);
            plot(taxispos, slice22, 'k-');
            xlim([min(taxispos) max(taxispos)]);
            ylim([min(slice22) max(slice22)]);
            xtick = [1 ceil(length(taxispos)/2) length(taxispos)];
            xtick = taxispos(xtick);
            xticklabel = xtick;
            xticklabel = ceil(xtick*10)/10;
            set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
            box on;
            tickpref;
            xlabel('Time (ms)');
            ylabel('Norm. Amplitude');
            title('STRF2 slice');


            subplot(5,1,5);
            hold on;
            plot(lags,c12,'k-');
            hb = bar(cc_delay,cc_q12);
            set(hb, 'facecolor', 'k');
            xtick = -50:5:50; %sort([-0 -25 -10 0 10 25 50]);
            set(gca,'xtick', xtick, 'xticklabel', xtick);
            xlim([min(lags) max(lags)]);
            xlim([-10 10]);
            ylim(1.03*[ min(min([c12(:); cc_q12(:)])) max(max([c12(:); cc_q12(:)]))]);
            box on;
            tickpref;
        %     legend('STRF Overlap', 'Spike Train Corr');
            xlabel('Delay (ms)');

            suptitle(sprintf('#%.0f: %s-site%.0f-chan%.0f model %.0f/%.0f %s', ...
                ii, exp, site, chan, model1(1), model2(1), stim));

            set(gcf,'position', [815 18 495 966]);

        end

    else
        fprintf('%s does not exist.\n', strfFile.name);
        ccoverlap(ii).overlap = [];

    end

    clear('spk', 'strf', 'trigger');

    pause;

end % (for ii)

return;











