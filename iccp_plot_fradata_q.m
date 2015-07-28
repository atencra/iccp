function iccpairs_plot_fradata_q(fradata)
% iccpairs_plot_fradata_bw ICC Pairwise FRA Q comparison
% 
%     iccpairs_plot_fradata_q(fradata)
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
Q = fradata.q;



for i = 1:4
   data = [Q(:,(i-1)*2+1) Q(:,i*2)];
   index = ~isnan(data(:,1)) & ~isnan(data(:,2));
   data = data(index,:);
   r = corrcoef(data(:,1), data(:,2));
   r = r(2);

   qdata{i} = data;
   qr{i} = r;
   qdiff{i} = abs(data(:,1) - data(:,2));
   qdiff_oct{i} = abs( log2(data(:,1)./data(:,2)) );

   qdiff_mn{i} = mean( qdiff{i} );
   qdiff_md{i} = median( qdiff{i} );
   qdiff_sd{i} = std( qdiff{i} );

   qdiff_oct_mn{i} = mean( qdiff_oct{i} );
   qdiff_oct_md{i} = median( qdiff_oct{i} );
   qdiff_oct_sd{i} = std( qdiff_oct{i} );
end % (for i)


minmin = 0.25;
maxmax = 8;

figure;
edges = linspace(0,4,21);

label = {'Q10', 'Q20', 'Q30', 'Q40'};
qtick = [0.25 0.5 1 2 4 8];

ydiffmax = [0.25 0.35 0.3 0.3];
ytick = {0:0.125:0.25, 0:0.175:1, 0:0.15:0.3, 0:0.15:0.3}; 
markersize = 1;

for i = 1:4

   data = qdata{i};
   datadiff = qdiff_oct{i};

   subplot(4,2,(i-1)*2+1);
   hold on
   xlim([minmin maxmax]);
   ylim([minmin maxmax]);
   plot(xlim, ylim, 'k-');
   [larger, smaller] = vs_largersmaller(qdata{i}(:,1),qdata{i}(:,2));
   scatter(larger,smaller,20, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])

   xlim([minmin maxmax])
   ylim([minmin maxmax])
   tickpref
   box off
   set(gca,'xtick', qtick, 'xticklabel', qtick);
   set(gca,'ytick', qtick, 'yticklabel', qtick);
   set(gca,'xscale', 'log', 'yscale', 'log');
   ylabel(sprintf('%s Neuron 2',label{i}))
   xlabel(sprintf('%s Neuron 1',label{i}))

   subplot(4,2,(i)*2);
   n = histc(qdiff_oct{i}, edges);
   n = n./sum(n);
   hb = bar(edges,n,'histc');
   set(hb, 'FaceColor', [0.6 0.6 0.6])
   xlabel('Q Difference (oct)')
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







