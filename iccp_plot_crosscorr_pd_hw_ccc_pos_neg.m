function iccp_plot_crosscorr_pd_hw_ccc_pos_neg(ccpairs, showstats)
% iccpairs_plot_crosscorr_pd_hw_ccc_pos_neg_hist Cross covariance parameters
% 
%    iccpairs_plot_crosscorr_pd_hw_ccc_pos_neg_hist(ccpairs, showstats)
% 
%    ccpairs : struct array, with each element holding the data for a pair
%    of ICC neurons.
% 
%    Splits data into 2 groups: 1) positive cross-covariance deflections
%    and 2) negative deflections.
% 
%    Plots peak delay, half-width, and cross correlation coefficient
%    for each data group.
%
%    If showstats = 1, summary statistics will be output. Otherwise, no
%    stats will be shown.

if ( nargin == 1 )
   showstats = 0;
end

close all;


data = iccp_get_ccpairs_pos_neg_data(ccpairs);

iccp_plot_ccpairs_pos_neg_pd_hw_ccc(data, showstats);

%iccp_plot_ccpairs_pos_pd_hw_ccc(data, showstats);

return;





function iccp_plot_ccpairs_pos_pd_hw_ccc(data, showstats)


sigfeature_pos = data.sigPos;

pd_pos = data.pdPos;
hw_pos = data.hwPos;
ccc_pos = data.cccPos;



edges_pd = linspace(0,5,11);
edges_hw = 0:0.5:6; %linspace(0,6,25);
edges_ccc = logspace(log10(0.001),log10(0.5),25);

count_pd_pos = histc(pd_pos(sigfeature_pos), edges_pd);
pdf_pd_pos = count_pd_pos / sum(count_pd_pos);

count_hw_pos = histc(hw_pos(sigfeature_pos), edges_hw);
pdf_hw_pos = count_hw_pos / sum(count_hw_pos);

count_ccc_pos = histc(ccc_pos(sigfeature_pos), edges_ccc);
pdf_ccc_pos = count_ccc_pos / sum(count_ccc_pos);


markersize = 4;

cmap = 0.5 * ones(1,3);

figure;

subplot(3,1,1)
hold on;
hp = plot(edges_pd, pdf_pd_pos, 'ko-', 'markersize', markersize, 'markerfacecolor', 'k', ...
'markeredgecolor', 'k');
set(hp, 'color', cmap(1,:));

xlabel('Peak Delay (ms)');
ylabel('Proportion')
xlim([min(edges_pd) max(edges_pd)]);
ylim([0 0.4]);
set(gca,'xtick', edges_pd(1:2:end), 'xticklabel', edges_pd(1:2:end));
ytick = 0:0.1:0.4;
set(gca,'ytick', ytick, 'yticklabel', ytick);
box off;
tickpref;
legend('Facilitative');
subplot_label(gca,'A');


subplot(3,1,2)
hold on;
hp = plot(edges_hw, pdf_hw_pos, 'ko-', 'markersize', markersize, 'markerfacecolor', 'k', ...
'markeredgecolor', 'k');
set(hp, 'color', cmap(1,:));

xlabel('Half Width (ms)')
ylabel('Proportion')
xlim([min(edges_hw) max(edges_hw)]);
ylim([0 0.5]);
set(gca,'ytick', 0:0.1:0.5, 'yticklabel', 0:0.1:0.5);
set(gca,'xtick', 0:6, 'xticklabel', 0:6);
box off;
tickpref;
subplot_label(gca,'B');



subplot(3,1,3)
hold on;
hp = plot(edges_ccc, pdf_ccc_pos, 'ko-', 'markersize', markersize, 'markerfacecolor', 'k', ...
'markeredgecolor', 'k');
set(hp, 'color', cmap(1,:));

set(gca,'xscale','log')
axis([0.000001 0.1 0 0.2])
ylim([0 0.21]);
ytick = 0:0.05:0.25;
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([min(edges_ccc) max(edges_ccc)]);
xlabel('CCC');
ylabel('Proportion');
box off
tickpref
subplot_label(gca,'C');

print_mfilename(mfilename);
set(gcf,'position', [1178 337 321 617]);


if ( showstats )
   f = simple_stats(pd_pos(sigfeature_pos));
   label = 'Facilitative';
   fprintf('%s Peak Delay:\n', label);
   fprintf('mean = %.3f\n', f.mn);
   fprintf('SD = %.3f\n', f.sd);
   fprintf('SE = %.3f\n', f.se);
   fprintf('median = %.3f\n', f.md);
   fprintf('MAD = %.3f\n', f.mad);
   fprintf('\n');

   f = simple_stats(hw_pos(sigfeature_pos));
   fprintf('\n');
   label = 'Facilitative';
   fprintf('%s Half-Width:\n', label);
   fprintf('mean = %.3f\n', f.mn);
   fprintf('SD = %.3f\n', f.sd);
   fprintf('SE = %.3f\n', f.se);
   fprintf('median = %.3f\n', f.md);
   fprintf('MAD = %.3f\n', f.mad);
   fprintf('\n');

   f = simple_stats(ccc_pos(sigfeature_pos));
   fprintf('\n');
   label = 'Facilitative';
   fprintf('%s CCC:\n', label);
   fprintf('mean = %.3f\n', f.mn);
   fprintf('SD = %.3f\n', f.sd);
   fprintf('SE = %.3f\n', f.se);
   fprintf('median = %.3f\n', f.md);
   fprintf('MAD = %.3f\n', f.mad);
   fprintf('\n');
end % (if)

return;




function iccp_plot_ccpairs_pos_neg_pd_hw_ccc(data, showstats)


sigfeature_pos = data.sigPos;

pd_pos = data.pdPos;
hw_pos = data.hwPos;
ccc_pos = data.cccPos;



sigfeature_neg = data.sigNeg;

pd_neg = data.pdNeg;
hw_neg = data.hwNeg;
ccc_neg = data.cccNeg;



edges_pd = linspace(0,5,11);
edges_hw = 0:0.5:6; %linspace(0,6,25);
edges_ccc = logspace(log10(0.001),log10(0.5),25);


count_pd_pos = histc(pd_pos(sigfeature_pos), edges_pd);
count_pd_neg = histc(pd_neg(sigfeature_neg), edges_pd);
pdf_pd_pos = count_pd_pos / sum(count_pd_pos);
pdf_pd_neg = count_pd_neg / sum(count_pd_neg);

fprintf('\n');
label = 'Peak delay';
iccp_ccc_pos_neg_summary(pd_pos(sigfeature_pos), pd_neg(sigfeature_neg), ...
   pd_pos(sigfeature_pos), label);


count_hw_pos = histc(hw_pos(sigfeature_pos), edges_hw);
count_hw_neg = histc(hw_neg(sigfeature_neg), edges_hw);
pdf_hw_pos = count_hw_pos / sum(count_hw_pos);
pdf_hw_neg = count_hw_neg / sum(count_hw_neg);

fprintf('\n');
label = 'Half Width';
iccp_ccc_pos_neg_summary(hw_pos(sigfeature_pos), hw_neg(sigfeature_neg), ...
   hw_pos(sigfeature_pos), label);


count_ccc_pos = histc(ccc_pos(sigfeature_pos), edges_ccc);
count_ccc_neg = histc(ccc_neg(sigfeature_neg), edges_ccc);
pdf_ccc_pos = count_ccc_pos / sum(count_ccc_pos);
pdf_ccc_neg = count_ccc_neg / sum(count_ccc_neg);

fprintf('\n');
label = 'CCC';
iccp_ccc_pos_neg_summary(ccc_pos(sigfeature_pos), ccc_neg(sigfeature_neg), ...
   ccc_pos(sigfeature_pos), label);


cmapcell = {'greens'};
markersize = 4;

for i = 1:length(cmapcell)

   cmap = brewmaps(cmapcell{i}, 4);
   cmap = cmap(1:3,:);

   %Population Histograms
   figure;
   subplot(3,1,1)
   hold on;
   hp = plot(edges_pd, pdf_pd_pos, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', cmap(1,:));
   set(hp, 'color', cmap(1,:));

   hp = plot(edges_pd, pdf_pd_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
    'markeredgecolor', cmap(3,:));
   set(hp, 'color', cmap(3,:));

   xlabel('Peak Delay (ms)');
   ylabel('Proportion')
   xlim([min(edges_pd) max(edges_pd)]);
   ylim([0 0.45]);
   set(gca,'xtick', edges_pd(1:2:end), 'xticklabel', edges_pd(1:2:end));
   ytick = 0:0.1:0.6;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   box off;
   tickpref;
   legend('Facilitative', 'Suppressive');
   subplot_label(gca,'A');


   subplot(3,1,2)
   hold on;
   hp = plot(edges_hw, pdf_hw_pos, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', cmap(1,:));
   set(hp, 'color', cmap(1,:));

   hp = plot(edges_hw, pdf_hw_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
    'markeredgecolor', cmap(3,:));
   set(hp, 'color', cmap(3,:));

   xlabel('Half Width (ms)')
   ylabel('Proportion')
   xlim([min(edges_hw) max(edges_hw)]);
   ylim([0 0.5]);
   set(gca,'ytick', 0:0.1:0.5, 'yticklabel', 0:0.1:0.5);
   set(gca,'xtick', 0:6, 'xticklabel', 0:6);
   box off;
   tickpref;
   subplot_label(gca,'B');



   subplot(3,1,3)
   hold on;
   hp = plot(edges_ccc, pdf_ccc_pos, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', cmap(1,:));
   set(hp, 'color', cmap(1,:));

   hp = plot(edges_ccc, pdf_ccc_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
    'markeredgecolor', cmap(3,:));
   set(hp, 'color', cmap(3,:));

   set(gca,'xscale','log')
   axis([0.000001 0.1 0 0.2])
   ylim([0 0.25]);
   ytick = 0:0.05:0.25;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   xlim([min(edges_ccc) max(edges_ccc)]);
   xlabel('CCC');
   ylabel('Proportion');
   box off
   tickpref
   subplot_label(gca,'C');


   print_mfilename(mfilename);
   set(gcf,'position', [1178 337 321 617]);

end %(for i)

if ( showstats )
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
end % (if)

return;




function data = iccp_get_ccpairs_pos_neg_data(ccpairs)


[pdPos, hwPos, cccPos, pdNeg, hwNeg, cccNeg, sigPos, sigNeg] = ...
    ccpairs_to_sigfeature(ccpairs);


data.sigPos = sigPos(:);
data.sigNeg = sigNeg(:);

data.pdPos = pdPos(:);
data.hwPos = hwPos(:);
data.cccPos = cccPos(:);

data.pdNeg = pdNeg(:);
data.hwNeg = hwNeg(:);
data.cccNeg = cccNeg(:);

return;
















