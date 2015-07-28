function [ccstr] = iccp_calc_spk_crosscorr(unitstr)
% iccp_calc_spk_crosscorr Cross-cov/corr from paired spike train data
% 
%     [ccpairs] = iccp_calc_spk_crosscorr(pairstrains);
% 
%     For the pairs of spike train data in pairstrains, estimate the 
%     correlation function for each pair of spike trains in pairstrains.
%
%     pairstrains is obtained from: [pairstrains] = iccpairs_get_icc_fra_spk_paired_data(spk);
% 
%     The data for each in pairstrains is returned in the struct array
%     ccpairs.
% 
%     [ccpairs] = iccp_calc_spk_crosscorr(pairstrains);

library('spike_train_cross_correlation');
library('agmon');

ccstr = unitstr;

for i = 1:length(ccstr)

    fprintf('Processing #%.0f of %.0f\n', i, length(ccstr) );

    fsspk = ccstr(i).fs;

    spiketimes1 = ccstr(i).spiketimes1;
    spiketimes2 = ccstr(i).spiketimes2;
    maxspiketime = max([ max(spiketimes1) max(spiketimes2) ] );
    spiketimes1 = [spiketimes1 maxspiketime];
    spiketimes2 = [spiketimes2 maxspiketime];

    spet1 = spiketimes1 ./ 1000 .* fsspk;
    spet2 = spiketimes2 ./ 1000 .* fsspk;

    fsd = 2000; % Hz = 0.5 ms resolution
    dt = 1 / fsd * 1000; % in ms

    train1 = spet2train(spet1, fsspk, fsd);
    train2 = spet2train(spet2, fsspk, fsd);


    % Compute cab = E[xa(t+u)xb(t)] = correlogram
    [r12, delay] = xcorr(train1, train2, 100); % make max lag = 50 * 0.5 ms = 25 ms
    r12 = real(r12);
    r12(find(r12 < 0)) = 0;
    delay = delay .* dt; % delay in ms, not bin number

    n1 = length(spet1);
    n2 = length(spet2);


    % estimate cross-covariance from correlogram
    [q12, conf_limit, stimdur] = cc_calc_cross_covariance(r12, n1, n2, dt);

    % peak delay, centroid, asymmetry, halfwidth, significance of
    % peak/correlogram
    [pd, cent, ca, hw, sigfeature] = ...
        cc_calc_cross_covariance_params(delay, q12, r12, conf_limit);


    % standard corr coefficient
    rho = cc_corrstrength(delay, r12, n1, n2, stimdur); % correlation coefficient

    % Agmon (2012) rate invariant cross-corr coeff
    cccval  = ccc(r12, spiketimes1, spiketimes2, dt);

    % Calculate binless correlogram
    a = spiketimes1 / 1000; % convert spiketimes to seconds
    b = spiketimes2 / 1000;
    binwidth = 1.5;
    [ccdelay, ccraw, ccflat, ccnorm] = cc_crosscorr_binless(a, b, binwidth);



    % save data to struct array
    ccstr(i).fsd = fsd;
    ccstr(i).dt = dt;
    ccstr(i).n1 = n1; % the number of spikes
    ccstr(i).n2 = n2; % the number of spikes
    ccstr(i).delay = delay; % delay for cross-cov/correlogram
    ccstr(i).r12 = r12; % correlogram
    ccstr(i).q12 = q12;
    ccstr(i).conf_limit = conf_limit;
    ccstr(i).rho = rho;
    ccstr(i).ccc = cccval;
    ccstr(i).peakdelay = pd;
    ccstr(i).centroid = cent;
    ccstr(i).asymmetry = ca;
    ccstr(i).halfwidth = hw;
    ccstr(i).significant = sigfeature;
    ccstr(i).ccdelay = ccdelay;
    ccstr(i).ccraw = ccraw;
    ccstr(i).ccflat = ccflat;
    ccstr(i).ccnorm = ccnorm;


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

    pause(0.1);

end % (for i)


return;


