function iccpairs_plot_mtfdata_btmf_bsmf(mtfdata)

% Pairs Data
rho = mtfdata.cc;
btmf = mtfdata.btmf;
bsmf = mtfdata.bsmf;

spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';
datatype = {'btmf', 'bsmf'};

xmin = [12.5 0.05];
xmax = [275 2.1];
ymin = xmin;
ymax = xmax;
xscale = {'log', 'log', 'log'};
yscale = {'log', 'log', 'log'};
edges = {linspace(0, 240, 17), linspace(0, 2, 21)};
tick = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
ticklabel = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
xlabel1 = {'bTMF (Hz)', 'bSMF (cyc/oct)'};
xlabel2 = {'bTMF Diff (Hz)', 'bSMF Diff (cyc/oct)'};

ymaxdiff = [.4 0.6];
ytick = {0:0.1:0.4, 0:0.2:0.6, linspace(0,0.6,4)}; 
xtick2 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1 2]}; 
ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1]}; 
cmap = brewmaps('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 3;
yax = {[0.004 1], [0.004 1]};
xax = {[0.1 400],[0.001 2]};

% x = {log10(logspace(log10(0.0001), log10(100), 1000)) };

x = {log10(logspace(log10(0.01), log10(100), 1000)), ...
log10(logspace(log10(0.001), log10(10), 1000))};

subplotnum = [1 3 5];

figure;

for i = 1:length(datatype)


    fignum = (i-1) + subplotnum;

   dataspktype = eval([datatype{i}]);
   [larger, smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
   datadiff = abs( dataspktype(:,1) - dataspktype(:,2) );

   subplot(3,2,fignum(1));
   hold on;
   plot(larger, smaller, 'ko', 'markerfacecolor', 0.6*ones(1,3), 'markersize', 2);
   plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');
   [r,p] = corrcoef(log10(dataspktype(:,1)), log10(dataspktype(:,2)));
   set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
   tickpref;
   set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
   set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
   xlim([xmin(i) xmax(i)]);
   ylim([ymin(i) ymax(i)]);
   xlabel(sprintf('%s %s', spktype_xlabel{1}, xlabel1{i}));
   ylabel(sprintf('%s %s', spktype_ylabel{1}, xlabel1{i}));
   title(sprintf('N = %.0f, r = %.3f, p = %.3f', size(dataspktype,1),...
      r(2), p(2)));



   subplot(3,2,fignum(2));
   n = histc(datadiff, edges{i});
   n = n ./ sum(n);
   hb = bar(edges{i}, n, 'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
   box off;
   tickpref; %set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   xtick = edges{i}(1:4:end);
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   range = max(edges{i}) - min(edges{i});
   xlim([min(edges{i}) max(edges{i})]);
   ylim([0 ymaxdiff(i)]);
   ylabel('Proportion');
   xlabel(xlabel2{i});



   subplot(3,2,fignum(3));
   hold on;
   plot(datadiff, rho, 'o', 'color', 'k', ...
   'markersize', markersize, 'markerfacecolor', 0.6*ones(1,3), ...
   'markeredgecolor', 'k');
   set(gca,'xtick', xtick2{i}, 'xticklabel', xtick2{i});
   set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
   tickpref;
   box off;
   xlim(xax{i});
   ylim(yax{i});
   xlabel(xlabel2{i});
   ylabel('Correlation Coefficient');
   [r,p] = corrcoef((datadiff(:)), log10(rho(:)));
   title(sprintf('r = %.3f, p = %.3f', r(1,2), p(1,2)));
   set(gca,'yscale', 'log');
   set(gca,'xscale', 'log');

end % for i

set(gcf,'position', [341 117 361 467]);
print_mfilename(mfilename);


return

