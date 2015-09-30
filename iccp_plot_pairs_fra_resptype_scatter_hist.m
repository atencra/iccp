function iccp_plot_pairs_fra_resptype_scatter_hist(nb_nb_ne, ne_nb_ne)
% iccp_plot_pairs_fra_resptype_scatter_hist Plot of response types for ICC pairs
% 
%     iccp_plot_pairs_fra_resptype_scatter_hist(nb_nb_ne, ne_nb_ne)
% 
%     Makes a scatter plot (x vs. y plot) for two response type measures. The
%     measures quantify the temporal profile of pure tone responses.
% 
%     nb_nb_ne : Nx2 array, each row represents one pair of neurons. Each column
%     represents the data for one neuron in a pair.
% 
%     The values are nb / (nb + ne), where 
%         nb = the number of spikes over the first half of tone stimulation, 
%     and ne = the number of spikes over the last half of tone stimulation.
% 
%     When nb/(nb+ne) ~ 1, then the neuron has a phasic/onset profile. When
%     nb/(nb+ne) ~ 0.5, the neuron has a tonic/sustained profile. There may be
%     rare occasions where nb/(nb+ne) ~ 0, which indicates a response when
%     the stimulus turns off.
% 
%     ne_nb_ne : Nx2 array, each row represents one pair of neurons. Each column
%     represents the data for one neuron in a pair. Now values near 0 indicate
%     onset neurons.



rtdata{1} = nb_nb_ne;
rtdata{2} = ne_nb_ne;

nb_nb_ne_diff = abs( nb_nb_ne(:,1) - nb_nb_ne(:,2) ); 
ne_nb_ne_diff = abs( ne_nb_ne(:,1) - ne_nb_ne(:,2) ); 

mn_nb_nb_ne_diff = mean( nb_nb_ne_diff );
md_nb_nb_ne_diff = median( nb_nb_ne_diff );
sd_nb_nb_ne_diff = std( nb_nb_ne_diff );

mn_ne_nb_ne_diff = mean( ne_nb_ne_diff );
md_ne_nb_ne_diff = median( ne_nb_ne_diff );
sd_ne_nb_ne_diff = std( ne_nb_ne_diff );

rtdiff{1} = nb_nb_ne_diff;
rtdiff{2} = ne_nb_ne_diff;

rtdiff_mn{1} = mn_nb_nb_ne_diff;
rtdiff_mn{2} = mn_ne_nb_ne_diff;

rtdiff_md{1} = md_nb_nb_ne_diff;
rtdiff_md{2} = md_ne_nb_ne_diff;

rtdiff_sd{1} = sd_nb_nb_ne_diff;
rtdiff_sd{2} = sd_ne_nb_ne_diff;


minmin = 0;
maxmax = 1;

edges = linspace(0,1,21);

label = {'NB-NB-NE', 'Phasic-Tonic Index'};
tick = [0:0.25:1];
nrows = length(label);

for i = 1:nrows

   figure;
   subplot(1,2,1);
   hold on
   xlim([minmin maxmax]);
   ylim([minmin maxmax]);
   plot(xlim, ylim, 'k-');
   [larger, smaller] = iccp_largersmaller(rtdata{i}(:,1),rtdata{i}(:,2));
   scatter(larger,smaller,10,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
   xlim([minmin maxmax])
   ylim([minmin maxmax])
   tickpref
   box off
   set(gca,'xtick', tick, 'xticklabel', tick);
   set(gca,'ytick', tick, 'yticklabel', tick);
   ylabel(sprintf('%s Neuron 2',label{i}))
   xlabel(sprintf('%s Neuron 1',label{i}))
   subplot_label(gca,'C');

   subplot(1,2,2);
   n = histc(rtdiff{i}, edges);
   n = n./sum(n);
   hb = bar(edges,n,'histc');
   set(hb, 'FaceColor', [0.6 0.6 0.6])
   xlabel('PTI Difference')
   xlim([min(edges) max(edges)]);
   ylim([0 0.4]);
   set(gca,'xtick', tick, 'xticklabel', tick);
   ytick = 0:0.1:0.7;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   tickpref
   box off
   subplot_label(gca,'D');


   [r,p] = corrcoef(larger, smaller );
   fprintf('\n');
   fprintf('%s\n',label{i});
   fprintf('r = %.4f, p = %.4f\n',r(2), p(2) );

   dstats = simple_stats(rtdiff{i});

   fprintf('\n');
   fprintf('%s Difference\n', label{i});
   fprintf('MN = %.4f, SD = %.4f, SE = %.4f, MD = %.4f, MAD = %.4f, n = %.0f\n',...
      dstats.mn, dstats.sd, dstats.se, dstats.md, dstats.mad, length(rtdata{i}(:,1)) );


   set(gcf,'position', [1255 657 438 161]);
   print_mfilename(mfilename);

end % (for i)


return;