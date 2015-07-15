function cc_plot_ccpairs(ccpairs)
% iccpairs_calc_spk_crosscorr Cross-cov/corr from spk data
% 
%     [ccpairs] = iccpairs_calc_spk_crosscorr(pairstrains);
% 
%     For the pairs of spike train data in pairstrains, estimate the 
%     correlation function for each spike train.
% 
%     The data for each in pairstrains is returned in the struct array
%     ccpairs.
% 
%     [ccpairs] = iccpairs_calc_spk_crosscorr(pairstrains);

library('spike_train_cross_correlation');
library('agmon');

close all;

for i = 1:length(ccpairs)

    fprintf('Processing #%.0f of %.0f\n', i, length(ccpairs) );

    fsd = ccpairs(i).fsd;
    dt = ccpairs(i).dt;
    n1 = ccpairs(i).n1;
    n2 = ccpairs(i).n2;
    delay = ccpairs(i).delay;
    r12 = ccpairs(i).r12;
    q12 = ccpairs(i).q12;
    conf_limit = ccpairs(i).conf_limit;
    rho = ccpairs(i).rho;
    cccval = ccpairs(i).ccc;
    pd = ccpairs(i).peakdelay;
    centroid = ccpairs(i).centroid;
    ca = ccpairs(i).asymmetry;
    halfwidth = ccpairs(i).halfwidth;
    significant = ccpairs(i).significant;
    ccdelay = ccpairs(i).ccdelay;
    ccraw = ccpairs(i).ccraw;
    ccflat = ccpairs(i).ccflat;
    ccnorm = ccpairs(i).ccnorm;


    % plot correlation functions as we go
    clf;
    subplot(3,1,1);
    bar(delay, r12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
    tickpref;
    box off;
    ylabel('Correlogram');


    % binless correlogram
    subplot(3,1,2);
    upper95qab = conf_limit;
    lower95qab = -conf_limit;
    hold on;
    ymin = min([min(q12) lower95qab]);
    ymax = max([max(q12) upper95qab]);
    range = ymax-ymin;
    bar(delay, q12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
    plot([0 0], [ymin-0.1*range ymax+0.1*range], 'k-');
    plot([pd pd], [ymin-0.1*range ymax+0.1*range], 'r-');
    plot([min(delay) max(delay)], [0 0], 'k-');
    plot([min(delay) max(delay)], [upper95qab upper95qab], 'r-');
    plot([min(delay) max(delay)], [lower95qab lower95qab], 'r-');
    tickpref;
    box off;
    ylabel('Cross-Covariance');


    % binless correlogram
    subplot(3,1,3);
    hold on;
    plot([0 0], [0 1.05*max(ccnorm)], '-', 'color', 0.6*ones(1,3), 'linewidth', 3);
    hb = bar(ccdelay, ccnorm);
    set(hb,'facecolor', 0*ones(1,3));
    set(hb,'edgecolor', 0*ones(1,3));
    set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
    xlim([-10 10]);
    xtick = ccdelay;
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    ylim([0 1.1*max(ccnorm)]);
    box off;
    xlabel('Delay (ms)');
    ylabel('CC Binless');

    set(gcf,'position', [569 94 560 892]);

    pause;

end % (for i)


return;


