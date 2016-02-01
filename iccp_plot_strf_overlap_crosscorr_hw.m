function iccp_plot_strf_overlap_crosscorr_hw(ccoverlap)
% iccp_plot_strf_overlap_crosscorr_hw  Compare STRF corr to spike train corr
%
% iccp_plot_strf_overlap_crosscorr_hw(ccoverlap)
% -----------------------------------------------------------
%
% ccoverlap : struct array holding the overlap data for each pair of neurons.
% ccoverlap has a field, named "overlap", that contains the calculations.
%
% ccoverlap is obtained from: [ccoverlap] = iccp_strf_overlap_crosscorr(ccpairs);
%


narginchk(1,1);

hw_overlap = [];
hw_cc = [];

for i = 1:length(ccoverlap)

    if ( ccoverlap(i).significant && ~isempty(ccoverlap(i).overlap) )

        hw_cc = [hw_cc ccoverlap(i).halfwidth];

        temp = ccoverlap(i).overlap.hw12;
        hw_overlap = [hw_overlap temp];

    end
end % (for i)


figure;
hold on;
plot(jitter(hw_cc,0.25), jitter(hw_overlap,0.25), ...
'ko', ...
'markersize', 4);
tickpref;
maxmax = max(max([hw_cc(:)' hw_overlap(:)']));
xlim([0 maxmax]);
ylim([0 maxmax]);
plot([0 maxmax], [0 maxmax], 'k-');
xlabel('Cross-Covariance Halfwidth (ms)');
ylabel('STRF Overlap Cross-Covariance Halfwidth (ms)');
title(mfilename);


[p,h] = ranksum(hw_cc, hw_overlap)


pdiff = 100 * abs(hw_cc - hw_overlap) ./ hw_cc;

simple_stats(pdiff,1);

return;







% Part 1.
subplot(4,3,1);
plot_strf_symmetric_colormap(strf1);
title('STRF1');

subplot(4,3,2);
plot_strf_symmetric_colormap(strf2);
title('STRF2');

subplot(4,3,4);
plot(taxis, slice11, 'k-');
xlim([min(taxis) max(taxis)]);
ylim([min(slice11) max(slice11)]);
box on;
tickpref;
title('STRF1 peak, STRF1 slice');

subplot(4,3,5);
plot(taxis, slice22, 'k-');
xlim([min(taxis) max(taxis)]);
ylim([min(slice22) max(slice22)]);
box on;
tickpref;
title('STRF2 peak, STRF2 slice');

subplot(4,3,6);
hold on;
patch([lag_low12 lag_low12 lag_high12 lag_high12],1*[-1 1 1 -1],0.85*[1 1 1]);
patch(1.25*[-1 -1 1 1],1*[-1 1 1 -1],[0.6 0.6 0.6]);
plot(lags,c12,'k-');
xtick = sort([-50 -25 -10 10 25 50]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
xlim([min(lags) max(lags)]);
box on;
tickpref;
title('Slice through each peak');


return;








