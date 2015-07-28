function iccpairs_plot_fradata_bw(fradata)
% iccpairs_plot_fradata_bw ICC Pairwise BW comparison
% 
%     iccpairs_plot_fradata_bw(fradata)
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
bw = fradata.bw;


for i = 1:4
   data = [bw(:,i) bw(:,i+4)];
   index = ~isnan(data(:,1)) & ~isnan(data(:,2));
   data = data(index,:);
   r = corrcoef(data(:,1), data(:,2));
   r = r(2);
   bwdata{i} = data;
   bwr{i} = r;
   bwdiff{i} = abs(data(:,1) - data(:,2));
   bwdiff_oct{i} = abs( log2(data(:,1)./data(:,2)) );
end % (for i)

% Plot of BW
figure
edges = linspace(0,4,21);
label = {'BW10', 'BW20', 'BW30', 'BW40'};
bwtick = [0.125 0.5 2 8 32];
ydiffmax = [0.25 0.35 0.35 0.35];
ytick = {0:0.125:0.25, 0:0.15:0.3, 0:0.15:0.3, 0:1.5:0.3}; 
minmin = 0.125;
maxmax = 32;

for i = 1:4
   data = bwdata{i};
   datadiff = bwdiff_oct{i};

   subplot(4,2,(i-1)*2+1);
   hold on
   xlim([minmin maxmax]);
   ylim([minmin maxmax]);
   plot(xlim, ylim, 'k-');
   [larger, smaller] = vs_largersmaller(bwdata{i}(:,1),bwdata{i}(:,2));
   scatter(larger,smaller,20,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
   xlim([minmin maxmax])
   ylim([minmin maxmax])
   tickpref
   box off
   set(gca,'xtick', bwtick, 'xticklabel', bwtick);
   set(gca,'ytick', bwtick, 'yticklabel', bwtick);
   set(gca,'xscale', 'log', 'yscale', 'log');
   ylabel(sprintf('%s Neuron 2 (kHz)',label{i}))
   xlabel(sprintf('%s Neuron 1 (kHz)',label{i}))


   subplot(4,2,(i)*2);
   n = histc(datadiff, edges);
   n = n./sum(n);
   hb = bar(edges,n,'histc');
   set(hb, 'FaceColor', [0.6 0.6 0.6])
   xlabel('BW Difference (oct)');
   ylabel('Proportion');
   xlim([min(edges) max(edges)]);
   ylim([0 ydiffmax(i)]);
   set(gca,'xtick', 0:1:4, 'xticklabel', 0:1:4);
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   tickpref
   box off


   fprintf('\n%s\n',label{i});
   [r,p] = corrcoef(larger, smaller);
   fprintf('r = %.4f, p = %.4f\n', r(2), p(2) );

   dstats = simple_stats(datadiff);

   fprintf('\n');
   fprintf('MN = %.4f, SD = %.4f, SE = %.4f, MD = %.4f, MAD = %.4f, n = %.0f\n',...
      dstats.mn, dstats.sd, dstats.se, dstats.md, dstats.mad, length(datadiff) );
   fprintf('\n');

end % (for i)

set(gcf,'position', [1003 248 361 551]);
print_mfilename(mfilename);


return;
















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
