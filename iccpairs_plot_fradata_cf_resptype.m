function iccpairs_plot_fradata_cf_resptype(fradata)
% iccpairs_plot_fradata_cf_resptype ICC Pairwise CF, resptype comparison
% 
%     iccpairs_plot_fradata_cf_resptype(fradata)
% 
% 
%     fradata : struct holding pairwise fra data. Obtained from
%     iccpairs_ccpairs_select_max_plot_fra_params.m
% 
%     Makes a 2x2 plot for CF and resptype. resptype is the phasic-tonic
%     index.
% 
%     Plot arranged as:
% 
%     A | B
%     --|--
%     C | D
% 
%     A : scatter plot of CF for each pair of neurons
%     B : histogram of CF differences
%     C : scatter plot of resptype for each pair of neurons
%     D : histogram of resptype differences


% Get data from struct
position = fradata.position;
cf = fradata.cf;
resptype = fradata.resptype;

% Calc differences and stats
cfdiff = abs( log2( cf(:,2)./cf(:,1) ) );
cfdstats = simple_stats(cfdiff);

resptype_diff = abs( resptype(:,1) - resptype(:,2) ); 
resptypedstats = simple_stats(resptype_diff);




% Plot CF1 vs CF2
cftick = [0.125 0.25 0.5 1 2 4 8 16 32];
figure;
subplot(2,2,1);
hold on
axis([0.2 32 0.2 32])
plot(xlim, ylim, 'k-');
[larger, smaller] = vs_largersmaller(cf(:,1), cf(:,2));
scatter(larger,smaller,20, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
tickpref
box off
set(gca,'xtick', cftick, 'xticklabel', cftick);
set(gca,'ytick', cftick, 'yticklabel', cftick);
set(gca,'xscale', 'log', 'yscale', 'log');
xlabel('CF (kHz) Neuron 1')
ylabel('CF (kHz) Neuron 2')
subplot_label(gca,'A');


% Plot CF differences histogram
subplot(2,2,2);
edges = linspace(0,1.5,13);
n = histc(cfdiff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
tickpref
box off
set(gca,'ytick', 0:0.25:0.5, 'yticklabel', 0:0.25:0.5);
set(gca,'xtick', 0:0.5:max(edges), 'xticklabel', 0:0.5:max(edges));
ylim([0 0.5])
xlim([min(edges) max(edges)]);
xlabel('CF Difference (oct)')
ylabel('Proportion');
subplot_label(gca,'B');


[r,p] = corrcoef(log10(larger), log10(smaller) );
fprintf('\n');
fprintf('CF Comparison\n');
fprintf('r = %.4f, p = %.4f\n',r(2), p(2) );

fprintf('\n');
fprintf('CF Difference (oct)\n');
fprintf('MN = %.4f, SD = %.4f, SE = %.4f, MD = %.4f, MAD = %.4f, n = %.0f\n',...
   cfdstats.mn, cfdstats.sd, cfdstats.se, cfdstats.md, cfdstats.mad, length(cfdiff) );



% Plot Resptype1 vs Resptype2
subplot(2,2,3);
hold on
minmin = 0;
maxmax = 1;
tick = [0:0.25:1];
xlim([minmin maxmax]);
ylim([minmin maxmax]);
plot(xlim, ylim, 'k-');
[larger, smaller] = vs_largersmaller(resptype(:,1),resptype(:,2));
scatter(larger,smaller,20,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
xlim([minmin maxmax])
ylim([minmin maxmax])
tickpref
box off
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
ylabel(sprintf('%s Neuron 2','PTI'))
xlabel(sprintf('%s Neuron 1','PTI'))
subplot_label(gca,'C');


% Plot Resptype1 vs Resptype2
subplot(2,2,4);
edges = linspace(0,0.6,13);
n = histc(resptype_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
xlim([min(edges) max(edges)]);
ylim([0 0.3]);
xtick = 0:0.2:0.6;
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = 0:0.1:0.3;
set(gca,'ytick', ytick, 'yticklabel', ytick);
tickpref
box off
xlabel('PTI Difference')
ylabel('Proportion');
subplot_label(gca,'D');


[r,p] = corrcoef(larger, smaller );
fprintf('\n');
fprintf('%s\n','PTI');
fprintf('r = %.4f, p = %.4f\n',r(2), p(2) );

fprintf('\n');
fprintf('%s Difference\n', 'PTI');
fprintf('MN = %.4f, SD = %.4f, SE = %.4f, MD = %.4f, MAD = %.4f, n = %.0f\n',...
resptypedstats.mn, resptypedstats.sd, resptypedstats.se, ...
resptypedstats.md, resptypedstats.mad, length(resptype_diff) );


set(gcf,'position', [1003 648 361 290]);
print_mfilename(mfilename);



return;
