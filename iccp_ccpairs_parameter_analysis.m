function iccp_ccpairs_parameter_analysis(ccpairs)
% ccpairs_parameter_analysis Analyzing various parameters from ccpairs data

% Looks at delay, half width, ccc, best frequency, q, latency, and
% bandwidth.

% Input:
% --------------------
% ccpairs: Struct array that is output from
% get_pairs_spkytype_strf_crosscorr. CURRENT DATA IS SAVED AS:
% icc_crosscorr_strf_new2

% Outputs:
% --------------------
% 1) Population histograms of peak delay, half width, and ccc, on both a
%       linear and a semilog scale.
% 2) Pairwise comparison scatterplots of best frequency, q, latency, and
%       bandwidth.
% 3) Comparison scatterplots of the differences in best frequency, q,
%       latency, and bandwidth, each against delay, half width, and ccc.


close all;


% Use new calculation of correlation coefficient if available
if ( isfield(ccpairs, 'ccc') )
   ccpairs_new = ccpairs;
else
   ccpairs_new = add_ccc_to_pairs_crosscorr_victor(ccpairs);
end

delay = abs([ccpairs_new.peakdelay]);
[min(delay) max(delay)];
width = [ccpairs_new.halfwidth];
ccc = [ccpairs_new.ccc];
index = find(ccc<0);
ccctemp = rand(size(index)) * 0.005;
ccc(index) = ccctemp;


BF1 = [ccpairs_new.cf1];
BF2 = [ccpairs_new.cf2];
Q1 = [ccpairs_new.q1];
Q2 = [ccpairs_new.q2];
L1 = [ccpairs_new.latency1];
L2 = [ccpairs_new.latency2];
BW1 = BF1./Q1;
BW2 = BF2./Q2;

deltaBF = abs(log2(BF1./BF2));
deltaQ = abs(log2(Q1./Q2));
deltaQ = abs(Q1 - Q2);
deltaL = abs(log2(L1./L2));
deltaL = abs(L1 - L2);
deltaBW = abs(log2(BW1./BW2));




iccpairs_plot_crosscorr_pd_hw_ccc_pos_neg_hist(ccpairs);


% vs_plot_delay_width_ccc_pophist(delay, width, ccc);

% vs_plot_bf_q_lat_bw_scatter(BF1, BF2, Q1, Q2, L1, L2, BW1, BW2);

% vs_plot_bf_q_lat_bw_scatter_diffhist(BF1, BF2, Q1, Q2, L1, L2, BW1, BW2);


% vs_plot_delta_pd_scatter(deltaBF, deltaQ, deltaL, deltaBW, width, delay);




% vs_plot_delta_width_scatter(deltaBF, deltaQ, deltaL, deltaBW, width);

% vs_plot_delta_ccc_scatter(deltaBF, deltaQ, deltaL, deltaBW, width, ccc);

% vs_plot_delta_pd_width_ccc_scatter(deltaBF, deltaQ, deltaL, deltaBW, delay, width, ccc);

return;





function iccpairs_plot_crosscorr_pd_hw_ccc_pos_neg_hist(ccpairs)

pd_pos = [ccpairs.pd_pos];
pd_pos = abs(pd_pos); % only absolute value matters, not direction
hw_pos = [ccpairs.hw_pos];
sigfeature_pos = [ccpairs.sigfeature_pos];
sigfeature_pos = logical(sigfeature_pos);
ccc_pos = [ccpairs.ccc_pos];
ccc_pos(ccc_pos < 0) = 0.0001;

pd_neg = [ccpairs.pd_neg];
pd_neg = abs(pd_neg);
hw_neg = [ccpairs.hw_neg];
sigfeature_neg = [ccpairs.sigfeature_neg];
sigfeature_neg = logical(sigfeature_neg);
ccc_neg = [ccpairs.ccc_neg];
ccc_neg(ccc_neg < 0) = 0.0001;


% Possibilities:
% Positive only peaks
% Negative only peaks
% Both Positive and Negative peaks
% Any significant peak


index_pos_only = sigfeature_pos & ~sigfeature_neg;
index_neg_only = ~sigfeature_pos & sigfeature_neg;
index_pos_neg = sigfeature_pos & sigfeature_neg;
index_any = sigfeature_pos | sigfeature_neg;


fprintf('Positive only peaks: %.0f\n', sum(index_pos_only) );
fprintf('Negative only peaks: %.0f\n', sum(index_neg_only) );
fprintf('Positive and Negative peaks: %.0f\n', sum(index_pos_neg) );
fprintf('Any peaks: %.0f\n', sum(index_any) );


edges_pd = 0:10; %linspace(0,25,11);
edges_hw = 0:0.5:6; %linspace(0,6,25);
edges_ccc = linspace(min(ccc_pos),max(ccc_pos),20);

e11 = 10.^(linspace(log(0.0001),log(max(pd_pos)),50));
e22 = 10.^(linspace(log(0.0001),log(max(hw_pos)),100));
% e33 = 10.^(linspace(log(0.0001),log(max(ccc)),50));
edges_ccc = logspace(log10(0.001),log10(0.5),25);


count_pd_pos = histc(pd_pos(sigfeature_pos), edges_pd);
count_pd_neg = histc(pd_neg(sigfeature_neg), edges_pd);
pdf_pd_pos = count_pd_pos / sum(count_pd_pos);
pdf_pd_neg = count_pd_neg / sum(count_pd_neg);

count_hw_pos = histc(hw_pos(sigfeature_pos), edges_hw);
count_hw_neg = histc(hw_neg(sigfeature_neg), edges_hw);
pdf_hw_pos = count_hw_pos / sum(count_hw_pos);
pdf_hw_neg = count_hw_neg / sum(count_hw_neg);

count_ccc_pos = histc(ccc_pos(sigfeature_pos), edges_ccc);
count_ccc_neg = histc(ccc_neg(sigfeature_neg), edges_ccc);
pdf_ccc_pos = count_ccc_pos / sum(count_ccc_pos);
pdf_ccc_neg = count_ccc_neg / sum(count_ccc_neg);


cmap = cschemes('gnbu', 4);
cmap = cmap(1:3,:);


%Population Histograms
figure;
subplot(3,1,1)
hold on;
hp = plot(edges_pd, pdf_pd_pos, 'o-', 'markersize', 3, 'markerfacecolor', cmap(1,:), ...
 'markeredgecolor', cmap(1,:));
set(hp, 'color', cmap(1,:));

hp = plot(edges_pd, pdf_pd_neg, 's-', 'markersize', 4, 'markerfacecolor', cmap(3,:), ...
 'markeredgecolor', cmap(3,:));
set(hp, 'color', cmap(3,:));
% hb1 = bar(edges_pd, pdf_pos, 'histc');
% set(hb1, 'FaceColor', [0.6 0.6 0.6]);
xlabel('Positive Peak Delay');
ylabel('Proportion')
title('Pos Peak Delay')
xlim([min(edges_pd) max(edges_pd)]);
ylim([0 0.6]);
set(gca,'xtick', edges_pd(1:2:end), 'xticklabel', edges_pd(1:2:end));
ytick = 0:0.2:0.6;
set(gca,'ytick', ytick, 'yticklabel', ytick);

box off;
tickpref;
legend('Facilitative', 'Suppressive');


subplot(3,1,2)
hold on;
hp = plot(edges_hw, pdf_hw_pos, 'o-', 'markersize', 3, 'markerfacecolor', cmap(1,:), ...
 'markeredgecolor', cmap(1,:));
set(hp, 'color', cmap(1,:));

hp = plot(edges_hw, pdf_hw_neg, 's-', 'markersize', 4, 'markerfacecolor', cmap(3,:), ...
 'markeredgecolor', cmap(3,:));
set(hp, 'color', cmap(3,:));
xlabel('Half Width')
ylabel('Proportion')
xlim([min(edges_hw) max(edges_hw)]);
ylim([0 0.5]);
set(gca,'ytick', 0:0.1:0.5, 'yticklabel', 0:0.1:0.5);
set(gca,'xtick', edges_hw, 'xticklabel', edges_hw);
box off;
tickpref;


subplot(3,1,3)
hold on;
hp = plot(edges_ccc, pdf_ccc_pos, 'o-', 'markersize', 3, 'markerfacecolor', cmap(1,:), ...
 'markeredgecolor', cmap(1,:));
set(hp, 'color', cmap(1,:));
hp = plot(edges_ccc, pdf_ccc_neg, 's-', 'markersize', 4, 'markerfacecolor', cmap(3,:), ...
 'markeredgecolor', cmap(3,:));
set(hp, 'color', cmap(3,:));
set(gca,'xscale','log')
axis([0.000001 0.1 0 0.2])
ylim([0 0.21]);
ytick = 0:0.05:0.25;
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([min(edges_ccc) max(edges_ccc)]);
xlabel('CCC');
ylabel('Proportion');
title('Pos CCC');
box off
tickpref

print_mfilename(mfilename);
set(gcf,'position', [155   154   321   617]);


f = simple_stats(pd_pos(sigfeature_pos));
s = simple_stats(pd_neg(sigfeature_neg));
label = 'Facilitative';
fprintf('%s Peak Delay:\n', label);
fprintf('mean = %.3f\n', f.mn);
fprintf('SD = %.3f\n', f.sd);
fprintf('SE = %.3f\n', f.se);
fprintf('median = %.3f\n', f.md);
fprintf('MAD = %.3f\n', f.mad);
fprintf('\n');
label = 'Suppressive';
fprintf('%s Peak Delay:\n', label);
fprintf('mean = %.3f\n', s.mn);
fprintf('SD = %.3f\n', s.sd);
fprintf('SE = %.3f\n', s.se);
fprintf('median = %.3f\n', s.md);
fprintf('MAD = %.3f\n', s.mad);


f = simple_stats(hw_pos(sigfeature_pos));
s = simple_stats(hw_neg(sigfeature_neg));
fprintf('\n');
label = 'Facilitative';
fprintf('%s Half-Width:\n', label);
fprintf('mean = %.3f\n', f.mn);
fprintf('SD = %.3f\n', f.sd);
fprintf('SE = %.3f\n', f.se);
fprintf('median = %.3f\n', f.md);
fprintf('MAD = %.3f\n', f.mad);
fprintf('\n');
label = 'Suppressive';
fprintf('%s Half-Width:\n', label);
fprintf('mean = %.3f\n', s.mn);
fprintf('SD = %.3f\n', s.sd);
fprintf('SE = %.3f\n', s.se);
fprintf('median = %.3f\n', s.md);
fprintf('MAD = %.3f\n', s.mad);


f = simple_stats(ccc_pos(sigfeature_pos));
s = simple_stats(ccc_neg(sigfeature_neg));
fprintf('\n');
label = 'Facilitative';
fprintf('%s CCC:\n', label);
fprintf('mean = %.3f\n', f.mn);
fprintf('SD = %.3f\n', f.sd);
fprintf('SE = %.3f\n', f.se);
fprintf('median = %.3f\n', f.md);
fprintf('MAD = %.3f\n', f.mad);
fprintf('\n');
label = 'Suppressive';
fprintf('%s CCC:\n', label);
fprintf('mean = %.3f\n', s.mn);
fprintf('SD = %.3f\n', s.sd);
fprintf('SE = %.3f\n', s.se);
fprintf('median = %.3f\n', s.md);
fprintf('MAD = %.3f\n', s.mad);


% fprintf('Half Width:\n');
% fprintf('mean = %.3f\n', nanmean(width));
% fprintf('SD = %.3f\n', nanstd(width));
% fprintf('median = %.3f\n', nanmedian(width));
% fprintf('MAD = %.3f\n', mad(ccc(~isnan(width) ) ));
% 
% fprintf('CCC:\n');
% fprintf('mean = %.3f\n', nanmean(ccc));
% fprintf('SD = %.3f\n', nanstd(ccc));
% fprintf('median = %.3f\n', nanmedian(ccc));
% fprintf('MAD = %.3f\n', mad(ccc(~isnan(ccc) ) ));

return;




function vs_plot_delay_width_ccc_pophist(delay, width, ccc)

% e1 = linspace(min(delay),max(delay),21);
edges_pd = 0:1.25:20; %linspace(0,25,11);
% e2 = linspace(min(width),max(width),20)
e2 = 0:0.5:6; %linspace(0,6,25);
e3 = linspace(min(ccc),max(ccc),20);

e11 = 10.^(linspace(log(0.0001),log(max(delay)),50));
e22 = 10.^(linspace(log(0.0001),log(max(width)),100));
% e33 = 10.^(linspace(log(0.0001),log(max(ccc)),50));
e33 = logspace(log10(0.001),log10(0.5),25);


x = histc(delay,edges_pd);
y = histc(width,e2);
z = histc(ccc,e3);

xx = histc(delay,e11);
yy = histc(width,e22);
zz = histc(ccc,e33);

x = x./sum(x);
y = y./sum(y);
z = z./sum(z);

xx = xx./sum(xx);
yy = yy./sum(yy);
zz = zz./sum(zz);

%Population Histograms
figure;
subplot(3,1,1)
hb1 = bar(edges_pd,x,'histc');
set(hb1, 'FaceColor', [0.6 0.6 0.6]);
xlabel('Peak Delay');
ylabel('Proportion')
title('Peak Delay Population Data')
xlim([min(edges_pd) max(edges_pd)]);
ylim([0 0.6]);
set(gca,'xtick', edges_pd(1:2:end), 'xticklabel', edges_pd(1:2:end));
box off;
tickpref;


subplot(3,1,2)
hb2 = bar(e2,y,'histc');
set(hb2, 'FaceColor', [0.6 0.6 0.6]);
xlabel('Half Width')
ylabel('Proportion')
xlim([min(e2) max(e2)]);
ylim([0 0.5]);
set(gca,'ytick', 0:0.1:0.5, 'yticklabel', 0:0.1:0.5);
set(gca,'xtick', e2, 'xticklabel', e2);
box off;
tickpref;


subplot(3,1,3)
hb33 = bar(e33,zz,'histc');
c11 = get(gca,'children');
set(c11(1),'marker','none');
set(gca,'xscale','log')
axis([0.000001 0.1 0 0.2])
set(hb33, 'FaceColor', [0.6 0.6 0.6]);
ylim([0 0.15]);
set(gca,'ytick', 0:0.05:0.15, 'yticklabel', 0:0.05:0.15);
xlim([min(e33) max(e33)]);
xlabel('CCC')
ylabel('Proportion')
box off
tickpref
print_mfilename(mfilename);
set(gcf,'position', [155   154   321   617]);


fprintf('Peak Delay:\n');
fprintf('mean = %.3f\n', nanmean(delay));
fprintf('SD = %.3f\n', nanstd(delay));
fprintf('median = %.3f\n', nanmedian(delay));
fprintf('MAD = %.3f\n', mad(ccc(~isnan(delay) ) ));

fprintf('Half Width:\n');
fprintf('mean = %.3f\n', nanmean(width));
fprintf('SD = %.3f\n', nanstd(width));
fprintf('median = %.3f\n', nanmedian(width));
fprintf('MAD = %.3f\n', mad(ccc(~isnan(width) ) ));

fprintf('CCC:\n');
fprintf('mean = %.3f\n', nanmean(ccc));
fprintf('SD = %.3f\n', nanstd(ccc));
fprintf('median = %.3f\n', nanmedian(ccc));
fprintf('MAD = %.3f\n', mad(ccc(~isnan(ccc) ) ));

return;



function vs_plot_bf_q_lat_bw_scatter(BF1, BF2, Q1, Q2, L1, L2, BW1, BW2) 


% Pairwise Plots
figure;

subplot(4,1,1);
hold on;
[larger, smaller] = vs_largersmaller(BF1,BF2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(BF1,BF2);
title(sprintf('BF1 vs BF2, r = %.2f' ,r(2)))
minmin = min([min(BF1) min(BF2)]);
maxmax = max([max(BF1) max(BF2)]);
minmin = 500;
maxmax = 20000;
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
set(gca,'xscale', 'log', 'yscale', 'log');
xlabel('Best Frequency (Hz) Neuron 1');
ylabel('Best Frequency (Hz) Neuron 2');
box off
tickpref;

subplot(4,1,2)
hold on;
[larger, smaller] = vs_largersmaller(Q1,Q2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(Q1,Q2);
title(sprintf('Q1 vs Q2, r = %.2f',r(2)))
minmin = min([min(Q1) min(Q2)]);
maxmax = max([max(Q1) max(Q2)]);
xlim([minmin maxmax]);
xlim([0.125 maxmax]);
ylim([minmin maxmax]);
ylim([0.125 maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
xtick = [0.125 0.25 0.5 1 2 4];
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'xscale', 'log', 'yscale', 'log');
xlabel('Spectral Tuning (Q) Neuron 1');
ylabel('Spectral Tuning (Q) Neuron 2');
box off
tickpref

subplot(4,1,3)
hold on;
[larger, smaller] = vs_largersmaller(L1,L2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(L1,L2);
title(sprintf('L1 vs L2, r = %.2f',r(2)))
minmin = min([min(L1) min(L2)]) - 0.5;
maxmax = max([max(L1) max(L2)]) + 0.5;
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
xlabel('Latency (ms) Neuron 1');
ylabel('Latency (ms) Neuron 2');
box off
tickpref

subplot(4,1,4)
hold on;
[larger, smaller] = vs_largersmaller(BW1,BW2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(BW1,BW2);
title(sprintf('BW1 vs BW2, r = %.2f',r(2)))
minmin = min([min(BW1) min(BW2)]);
maxmax = max([max(BW1) max(BW2)]);
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
xlabel('Bandwidth (Hz) Neuron 1');
ylabel('Bandwidth (Hz) Neuron 2');
set(gca,'xscale', 'log', 'yscale', 'log');
box off
tickpref
print_mfilename(mfilename);
set(gcf,'position', [141   101   321   880]);

return;



function vs_plot_bf_q_lat_bw_scatter_diffhist(BF1, BF2, Q1, Q2, L1, L2, BW1, BW2) 

deltaBF = abs(log2(BF1./BF2));
deltaQ = abs(log2(Q1./Q2));
deltaQ = abs(Q1 - Q2);
deltaL = abs(log2(L1./L2));
deltaL = abs(L1 - L2);
deltaBW = abs(log2(BW1./BW2));



% Pairwise Plots
figure;

subplot(4,2,1);
hold on;
[larger, smaller] = vs_largersmaller(BF1,BF2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(BF1,BF2);
title(sprintf('BF1 vs BF2, r = %.2f' ,r(2)))
minmin = min([min(BF1) min(BF2)]);
maxmax = max([max(BF1) max(BF2)]);
minmin = 500;
maxmax = 20000;
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
set(gca,'xscale', 'log', 'yscale', 'log');
xlabel('Best Frequency (Hz) Neuron 1');
ylabel('Best Frequency (Hz) Neuron 2');
box off
tickpref;

subplot(4,2,2)
edges = 0:0.1:2;
x = histc(deltaBF,edges);
x = x ./ sum(x);
hb = bar(edges,x,'histc');
children = get(gca,'children');
set(children(1),'marker','none');
set(hb, 'FaceColor', [0.6 0.6 0.6]);
xtick = 0:0.5:2;
ytick = 0:0.1:0.6;
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
ylim([0 0.6]);
xlim([min(edges) max(edges)]);
xlabel('BF Difference (oct)')
ylabel('Proportion')
box off
tickpref



subplot(4,2,3)
hold on;
[larger, smaller] = vs_largersmaller(Q1,Q2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(Q1,Q2);
title(sprintf('Q1 vs Q2, r = %.2f',r(2)))
minmin = min([min(Q1) min(Q2)]);
maxmax = max([max(Q1) max(Q2)]);
xlim([minmin maxmax]);
xlim([0.125 maxmax]);
ylim([minmin maxmax]);
ylim([0.125 maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
xtick = [0.125 0.25 0.5 1 2 4];
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'xscale', 'log', 'yscale', 'log');
xlabel('Spectral Tuning (Q) Neuron 1');
ylabel('Spectral Tuning (Q) Neuron 2');
box off
tickpref


subplot(4,2,4)
edges = 0:0.2:4;
x = histc(deltaQ,edges);
x = x ./ sum(x);
hb = bar(edges,x,'histc');
children = get(gca,'children');
set(children(1),'marker','none');
set(hb, 'FaceColor', [0.6 0.6 0.6]);
xtick = 0:1:4;
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = 0:0.1:0.3;
set(gca,'ytick', ytick, 'yticklabel', ytick);
ylim([0 0.3]);
xlim([min(edges) max(edges)]);
xlabel('Q Difference')
ylabel('Proportion')
box off
tickpref


subplot(4,2,5)
hold on;
[larger, smaller] = vs_largersmaller(L1,L2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(L1,L2);
title(sprintf('L1 vs L2, r = %.2f',r(2)))
minmin = min([min(L1) min(L2)]) - 0.5;
maxmax = max([max(L1) max(L2)]) + 0.5;
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
xlabel('Latency (ms) Neuron 1');
ylabel('Latency (ms) Neuron 2');
box off
tickpref


subplot(4,2,6)
edges = 0:0.5:10;
x = histc(deltaL,edges);
x = x ./ sum(x);
hb = bar(edges,x,'histc');
children = get(gca,'children');
set(children(1),'marker','none');
set(hb, 'FaceColor', [0.6 0.6 0.6]);
xtick = 0:2:10;
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = 0:0.1:0.3;
set(gca,'ytick', ytick, 'yticklabel', ytick);
ylim([0 0.31]);
xlim([min(edges) max(edges)]);
xlabel('Latency Difference(ms)')
ylabel('Proportion')
box off
tickpref


subplot(4,2,7)
hold on;
[larger, smaller] = vs_largersmaller(BW1,BW2);
plot(larger,smaller,'ko', 'markerfacecolor', 'k', 'markersize', 2);
r = corrcoef(BW1,BW2);
title(sprintf('BW1 vs BW2, r = %.2f',r(2)))
minmin = min([min(BW1) min(BW2)]);
maxmax = max([max(BW1) max(BW2)]);
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot([minmin maxmax], [minmin maxmax], 'k-');
xlabel('Bandwidth (Hz) Neuron 1');
ylabel('Bandwidth (Hz) Neuron 2');
set(gca,'xscale', 'log', 'yscale', 'log');
box off
tickpref
print_mfilename(mfilename);
set(gcf,'position', [141   101   321   880]);

subplot(4,2,8)
edges = 0:0.2:5;
x = histc(deltaBW,edges);
x = x ./ sum(x);
hb = bar(edges,x,'histc');
children = get(gca,'children');
set(children(1),'marker','none');
set(hb, 'FaceColor', [0.6 0.6 0.6]);
xtick = 0:1:5;
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = 0:0.1:0.3;
set(gca,'ytick', ytick, 'yticklabel', ytick);
ylim([0 0.3]);
xlim([min(edges) max(edges)]);
xlabel('BW Difference (oct)')
ylabel('Proportion')
box off
tickpref
print_mfilename(mfilename);
set(gcf,'position', [141   231   438   750]);

return;



function vs_plot_delta_pd_scatter(deltaBF, deltaQ, deltaL, deltaBW, width, delay)

index = ~isnan(width); % Only process if we could calculate Half-width

% Parameter Comparisons
figure;
subplot(4,1,1)
plot(deltaBF(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta BF vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DBF (oct)');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(4,1,2)
plot(deltaQ(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta Q vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DQ');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(4,1,3)
plot(deltaL(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta L vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DL (ms)');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(4,1,4)
plot(deltaBW(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta BW vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DBW (oct)');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;
set(gcf,'position', [137 52 321 918]);

return;



function vs_plot_delta_width_scatter(deltaBF, deltaQ, deltaL, deltaBW, width)

index = ~isnan(width); % Only process if we could calculate Half-width

figure;
subplot(4,1,1)
plot(deltaBF(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta BF vs Half Width')
xlabel('DBF (oct)');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(4,1,2)
plot(deltaQ(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta Q vs Half Width')
xlabel('DQ');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(4,1,3)
plot(jitter(deltaL(index)),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta L vs Half Width')
xlabel('DL (ms)');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(4,1,4)
plot(deltaBW(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta L vs Half Width')
xlabel('DBW (oct)');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;
set(gcf,'position', [137 52 321 918]);

return;



function vs_plot_delta_ccc_scatter(deltaBF, deltaQ, deltaL, deltaBW, width, ccc)

index = ~isnan(width); % Only process if we could calculate Half-width

figure;
subplot(4,1,1)
plot(deltaBF(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta BF vs CCC')
xlabel('DBF (oct)');
ylabel('CCC');
ylim([0.001 1]);
box off
tickpref

subplot(4,1,2)
plot(deltaQ(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta Q vs CCC')
xlabel('DQ');
ylabel('CCC');
ylim([0.001 1]);
box off
tickpref

subplot(4,1,3)
plot(jitter(deltaL(index)),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta L vs CCC')
xlabel('DL (ms)');
ylabel('CCC');
ylim([0.001 1]);
box off
tickpref

subplot(4,1,4)
plot(deltaBW(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta BW vs CCC')
xlabel('DBW (oct)');
ylabel('CCC');
ylim([0.001 1]);
box off
tickpref
print_mfilename(mfilename);
set(gcf,'position', [137 52 321 918]);



return;



function vs_plot_delta_pd_width_ccc_scatter(deltaBF, deltaQ, deltaL, deltaBW, delay, width, ccc)


% [deltaBF(:) delay(:) width(:) ccc(:)]

index = ~isnan(width); % Only process if we could calculate Half-width

% Parameter Comparisons
figure;
subplot(3,4,1)
plot(deltaBF(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta BF vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DBF (oct)');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(3,4,5)
plot(deltaBF(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta BF vs Half Width')
xlabel('DBF (oct)');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(3,4,9)
plot(deltaBF(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta BF vs CCC')
xlabel('DBF (oct)');
ylabel('CCC');
ylim([0.0001 1]);
box off
tickpref


subplot(3,4,2)
plot(deltaQ(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta Q vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DQ');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(3,4,6)
plot(deltaQ(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta Q vs Half Width')
xlabel('DQ');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(3,4,10)
plot(deltaQ(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta Q vs CCC')
xlabel('DQ');
ylabel('CCC');
ylim([0.0001 1]);
box off
tickpref


subplot(3,4,3)
plot(deltaL(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta L vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DL (ms)');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(3,4,7)
plot(deltaL(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta L vs Half Width')
xlabel('DL (ms)');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(3,4,11)
plot(deltaL(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta L vs CCC')
xlabel('DL (ms)');
ylabel('CCC');
ylim([0.0001 1]);
box off
tickpref


subplot(3,4,4)
plot(deltaBW(index), delay(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title(sprintf('Delta L vs Peak Delay; n = %.0f', length(delay(index))));
xlabel('DBW (oct)');
ylabel('Peak Delay (ms)');
ylim([0 20]);
tickpref;
box off;

subplot(3,4,8)
plot(deltaBW(index),width(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
title('Delta L vs Half Width')
xlabel('DBW (oct)');
ylabel('Half Width (ms)');
ylim([0 6]);
tickpref;
box off;

subplot(3,4,12)
plot(deltaBW(index),ccc(index),'ko', 'markerfacecolor', 'k', 'markersize', 2);
set(gca,'yscale', 'log');
% set(gca,'xscale', 'log');
title('Delta BW vs CCC')
xlabel('DBW (oct)');
ylabel('CCC');
ylim([0.0001 1]);
box off
tickpref
print_mfilename(mfilename);

return;


