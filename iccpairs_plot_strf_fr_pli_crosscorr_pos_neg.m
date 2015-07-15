function iccpairs_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs)
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% -----------------------------------------------------------------------
%
%
%
% caa 4/21/11

close all;

pd_pos = [ccpairs.pd_pos];
hw_pos = [ccpairs.hw_pos];
sigfeature_pos = [ccpairs.sigfeature_pos];
sigfeature_pos = logical(sigfeature_pos);
ccc_pos = [ccpairs.ccc_pos];
ccc_pos(ccc_pos < 0) = 0.0001;

pd_neg = [ccpairs.pd_neg];
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


% Data
fr1 = [ccpairs.fr1];
fr2 = [ccpairs.fr2];

pli1 = [ccpairs.pli1];
pli2 = [ccpairs.pli2];

si1 = [ccpairs.sepindex1];
si2 = [ccpairs.sepindex2];

fr = [fr1(:) fr2(:)];
pli = [pli1(:) pli2(:)];
si = [si1(:) si2(:)];

ccpairs_fr = fr;
ccpairs_pli = pli;
ccpairs_si = si;

xmin = [0.1 0.01 0.2];
xmax = [30 1 1];
ymin = xmin;
ymax = xmax;
xscale = {'log', 'log', 'log'};
yscale = {'log', 'log', 'log'};
edges = {0:2.5:20, 0:0.05:0.4, 0:0.075:0.6};
edges = {linspace(0,100,21), linspace(0,0.4,17), 0:0.075:0.6};
tick = {[0.1 1 10 30], [0.01 0.1 1], [0:0.2:1]};
ticklabel = {[5 10 20 30], [0.5 1 2 4 8], [0:10:50]};
xlabel1 = {'Firing Rate', 'RPI', 'Separability Index'};
xlabel2 = {'FR Difference (Hz)', 'RPI Difference', 'Separability Index Diff'};


spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'fr', 'pli', 'si'};
celltype = {'ccpairs'};

dataxlabel = {'Firing Rate (Hz) Neuron 1', 'RPI Neuron 1' };
dataylabel = {'Firing Rate (Hz) Neuron 2', 'RPI Neuron 2' };

datadiffxlabel = {'FR Difference (Hz)','RPI Difference (Hz)'};
datadiffylabel = {'Cross Correlation Coefficient','Cross Correlation Coefficient' };

ymaxdiff = [0.75 0.6 0.65];
ytick = {0:0.1:0.7, linspace(0,0.6,5), linspace(0,0.6,4)}; 
xtick2 = {0:20:100, 0:0.1:0.4, [0.0001 0.001 0.01 0.1 1]}; 
xtick3 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 4;
yax = {[0.001 1], [0.001 1], [0.0001 2]};
xax = {[0.01 100],[0.0001 1] [0.0001 2]};

% x = {log10(logspace(log10(0.0001), log10(100), 1000)) };

x = {log10(logspace(log10(0.01), log10(100), 1000)), ...
log10(logspace(log10(0.0001), log10(1), 1000))};

subplotnum = [1 3 5];

for i = 1:2 %length(datatype)


   figure;
   dataspktype = eval(['ccpairs_' datatype{i}]);
   dataspktype_pos = dataspktype(index_pos_only,:);
   dataspktype_neg = dataspktype(index_neg_only,:);
   dataspktype_pos_neg = dataspktype(index_pos_neg,:);

   [larger,smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
   datadiff = abs( larger - smaller );

% [length(datadiff) length(ccc_pos) length(ccc_neg)]

   [larger_pos,smaller_pos] = vs_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));
   datadiff_pos = abs( larger_pos - smaller_pos );

   [larger_neg,smaller_neg] = vs_largersmaller(dataspktype_neg(:,1),dataspktype_neg(:,2));
   datadiff_neg = abs( larger_neg - smaller_neg );

   [larger_pos_neg,smaller_pos_neg] = vs_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));
   datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );

   iccpairs_ccc_pos_neg_summary(datadiff_pos, datadiff_neg, datadiff_pos_neg, datatype{i});


   subplot(3,1,1);
   hold on;
   plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');

   cmap = cschemes('blues', 4);
   cmap = cmap(1:3,:);

   plot(larger_pos, smaller_pos, 'o', 'color', cmap(1,:), ...
      'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
      'markeredgecolor', cmap(1,:));

   plot(larger_neg, smaller_neg, 'o', 'color', cmap(2,:), ...
      'markersize', markersize, 'markerfacecolor', cmap(2,:), ...
      'markeredgecolor', cmap(2,:));

   plot(larger_pos_neg, smaller_pos_neg, 'o', 'color', cmap(3,:), ...
      'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
      'markeredgecolor', cmap(3,:));

   set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
   tickpref;
   set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
   set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
   xlim([xmin(i) xmax(i)]);
   ylim([ymin(i) ymax(i)]);
%    xlabel(sprintf('%s %s', xlabel1{i}, spktype_xlabel{1}));
%    ylabel(sprintf('%s %s', xlabel1{i}, spktype_ylabel{1}));
   xlabel( dataxlabel{i} );
   ylabel( dataylabel{i} );
   subplot_label(gca,'A');



   subplot(3,1,2);
   hold on;
   n = histc(datadiff_pos, edges{i});
   pdf_pos = n ./ sum(n);
   hp = plot(edges{i}, pdf_pos, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
   set(hp, 'color', cmap(1,:));

   n = histc(datadiff_neg, edges{i});
   pdf_neg = n ./ sum(n);
   hp = plot(edges{i}, pdf_neg, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
   set(hp, 'color', cmap(2,:));

   n = histc(datadiff_pos_neg, edges{i});
   pdf_pos_neg = n ./ sum(n);
   hp = plot(edges{i}, pdf_pos_neg, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
   set(hp, 'color', cmap(3,:));

   hold on;
   box off;
   tickpref;
   xtick = edges{i}(1:2:end);
   set(gca,'xtick', xtick2{i}, 'xticklabel', xtick2{i});
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   range = max(edges{i}) - min(edges{i});
   xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
   ylim([0 ymaxdiff(i)]);
   ylabel('Proportion');
   xlabel(xlabel2{i});
   legend('Pos', 'Neg', 'Pos+Neg');
   subplot_label(gca,'B');


     
   subplot(3,1,3);
   hold on;

   cmap = cschemes('greens', 4);
   cmap = cmap(1:3,:);

   plot(datadiff(sigfeature_pos), ccc_pos(sigfeature_pos), 'o', 'color', cmap(1,:), ...
      'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
      'markeredgecolor', cmap(1,:));

   [r,p] = corrcoef( log10(datadiff(sigfeature_pos)), log10(ccc_pos(sigfeature_pos)) );
   fprintf('%s: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));


   plot(datadiff(sigfeature_neg), ccc_neg(sigfeature_neg), 'o', 'color', cmap(3,:), ...
      'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
      'markeredgecolor', cmap(3,:));

   [r,p] = corrcoef( log10(datadiff(sigfeature_neg)), log10(ccc_neg(sigfeature_neg)) );
   fprintf('%s: Sup: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));


   legend('Facilitative', 'Suppressive');

   set(gca,'xtick', xtick3{i}, 'xticklabel', xtick3{i});
   set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
   tickpref;
   box off;
   xlim(xax{i});
   ylim(yax{i});
   xlabel(xlabel2{i});
   ylabel('Correlation Coefficient');
   set(gca,'yscale', 'log');
   set(gca,'xscale', 'log');
   subplot_label(gca,'C');


   fprintf('\n');
   fprintf('%s\n', datatype{i});
   [r,p] = corrcoef( log10(larger), log10(smaller) );
   fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
   fprintf('%s diff\n', datatype{i});
   x = datadiff(sigfeature_pos);
   fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
      length(x), mean(x), std(x), std(x)./sqrt(length(x)));
   fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
      median(datadiff), mad(datadiff), max(datadiff));
   x = datadiff(sigfeature_pos);
   y = ccc_pos(sigfeature_pos);
   [r,p] = corrcoef(x(:), log10(y(:)));
   fprintf('ccc vs diff, r = %.4f, p = %.4f', r(2), p(2));
   fprintf('\n');

   set(gcf,'position', [1434 131 321 694]);
   print_mfilename(mfilename);

end % for i




return;

