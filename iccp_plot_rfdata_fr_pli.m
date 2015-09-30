function iccp_plot_rfdata_fr_pli(rfdata)
% iccp_plot_rfdata_fr_pli(rfdata) Plot firing rate and phase-locking index
%


rho = rfdata.cc;

fr = rfdata.fr;
fr1 = fr(:,1);
fr2 = fr(:,2);

pli = rfdata.pli;
pli1 = pli(:,1);
pli2 = pli(:,2);


% ccpairs_rho = [ccpairs.rho];
% ccpairs_rho = [ccpairs.ccc];
% ccpairs_rho(ccpairs_rho<=0) = 0.0001;
% ccpairs_sig = logical([ccpairs.significant]);


% ccpairs_fr = fr;
% ccpairs_pli = pli;
% ccpairs_si = si;


xmin = [0.1 0.01 0.2];
xmax = [100 1 1];
ymin = xmin;
ymax = xmax;
xscale = {'log', 'log', 'log'};
yscale = {'log', 'log', 'log'};
edges = {0:2.5:20, 0:0.05:0.4, 0:0.075:0.6};
edges = {linspace(0,80,17), linspace(0,0.4,17), 0:0.075:0.6};
tick = {[0.1 1 10 100], [0.01 0.1 1], [0:0.2:1]};
ticklabel = {[5 10 20 30], [0.5 1 2 4 8], [0:10:50]};
xlabel1 = {'Firing Rate', 'RPI', 'Separability Index'};
xlabel2 = {'FR Difference (Hz)', 'RPI Difference', 'Separability Index Diff'};


spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'fr', 'pli'};
celltype = {'ccpairs'};

dataxlabel = {'Firing Rate (Hz) Neuron 1', 'RPI Neuron 1' };
dataylabel = {'Firing Rate (Hz) Neuron 2', 'RPI Neuron 2' };

datadiffxlabel = {'FR Difference (Hz)','RPI Difference (Hz)'};
datadiffylabel = {'Cross Correlation Coefficient','Cross Correlation Coefficient' };

ymaxdiff = [0.85 0.41 0.65];
ytick = {0:0.2:0.8, linspace(0,0.4,5), linspace(0,0.6,4)}; 
xtick2 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
cmap = brewmaps('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 3;
yax = {[0.001 1], [0.001 1], [0.0001 2]};
xax = {[0.01 100],[0.0001 1] [0.0001 2]};

% x = {log10(logspace(log10(0.0001), log10(100), 1000)) };

x = {log10(logspace(log10(0.01), log10(100), 1000)), ...
log10(logspace(log10(0.0001), log10(1), 1000))};

subplotnum = [1 3 5];

figure;

for i = 1:2 %length(datatype)

    fignum = (i-1) + subplotnum;

    dataspktype = eval([datatype{i}]);
    %    rho = eval(['ccpairs_rho']);
    %    sig = eval(['ccpairs_sig']);
    [larger,smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
    datadiff = abs( larger - smaller );

    % Scatter plot of values
    subplot(3,2,fignum(1));
    hold on;
    plot(larger, smaller, ...
    'ko', 'markerfacecolor', 'k', 'markersize', markersize);
    plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');
    set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
    tickpref;
    set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
    set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
    xlim([xmin(i) xmax(i)]);
    ylim([ymin(i) ymax(i)]);
    xlabel(sprintf('%s %s', xlabel1{i}, spktype_xlabel{1}));
    ylabel(sprintf('%s %s', xlabel1{i}, spktype_ylabel{1}));
    x = xlim;
    y = ylim;
    text(x(1), 1.05*y(2), 'A');


    % Difference between values for each pair of neurons
    subplot(3,2,fignum(2));
    n = histc(datadiff, edges{i});
    n = n ./ sum(n);
    hb = bar(edges{i},n,'histc');
    set(hb, 'facecolor', [0.7 0.7 0.7]);
    hold on;
    box off;
    tickpref;
    xtick = edges{i}(1:4:end);
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
    range = max(edges{i}) - min(edges{i});
    xlim([min(edges{i}) max(edges{i})]);
    ylim([0 ymaxdiff(i)]);
    ylabel('Proportion');
    xlabel(xlabel2{i});
    x = xlim;
    y = ylim;
    text(x(1), 1.05*y(2), 'B');


    % Correlation coefficient vs. difference of firing rate
    subplot(3,2,fignum(3));
    hold on;
    plot(datadiff, rho, 'o', 'color', 'k', ...
    'markersize', markersize, 'markerfacecolor', 'k', ...
    'markeredgecolor', 'k');
    set(gca,'xtick', xtick2{i}, 'xticklabel', xtick2{i});
    set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
    tickpref;
    box off;
    xlim(xax{i});
    ylim(yax{i});
    xlabel(xlabel2{i});
    ylabel('Correlation Coefficient');
    set(gca,'yscale', 'log');
    set(gca,'xscale', 'log');
    x = xlim;
    y = ylim;
    text(x(1), 1.05*y(2), 'C');


    fprintf('\n');
    fprintf('%s\n', datatype{i});
    [r,p] = corrcoef( log10(larger), log10(smaller) );
    fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
    fprintf('%s diff\n', datatype{i});
    fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
    length(datadiff),mean(datadiff), std(datadiff),...
    std(datadiff)./sqrt(length(datadiff)));
    fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
    median(datadiff), mad(datadiff), max(datadiff));
    [r,p] = corrcoef((datadiff(:)), log10(rho(:)));
    fprintf('ccc vs diff, r = %.4f, p = %.4f', r(2), p(2));
    fprintf('\n');

end % for i

set(gcf,'position',[247   101   361 551]);
print_mfilename(mfilename);

% se = sd ./ sqrt(length(datadiff));

return;










