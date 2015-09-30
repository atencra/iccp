function iccp_plot_pairs_fra_resptype_diff_rand(nb_nb_ne, ne_nb_ne)
% iccp_plot_pairs_fra_resptype_diff_rand Histogram of response type differences
% 
%     iccp_plot_pairs_fra_resptype_diff_rand(nb_nb_ne, ne_nb_ne)
%
%     Makes a histogram of response type differences for pairs of neurons.
%     Also performs a randomization test to determine if the distribution
%     of differences could be expected by a random ordering of distribution
%     of neurons.
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

figure;
edges = linspace(0,1,21);

label = {'NB-NB-NE', 'NE-NB-NE'};
tick = [0:0.25:1];
nrows = length(label);

for i = 1:nrows

   [larger, smaller] = iccp_largersmaller(rtdata{i}(:,1),rtdata{i}(:,2));

   subplot(nrows,2,(i-1)*2+1);
   n = histc(rtdiff{i}, edges);
   n = n./sum(n);
   hb = bar(edges,n,'histc');
   set(hb, 'FaceColor', [0.6 0.6 0.6])
   xlabel('Response Type Difference')
   ylabel('Proportion');
   xlim([min(edges) max(edges)]);
   ylim([0 0.4]);
   set(gca,'xtick', tick, 'xticklabel', tick);
   ytick = 0:0.1:0.7;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   tickpref
   box off

   [r,p] = corrcoef(larger, smaller );
   fprintf('\n');
   fprintf('%s\n',label{i});
   fprintf('r = %.4f, p = %.4f\n',r(2), p(2) );

   fprintf('\n');
   fprintf('%s Difference\n', label{i});
   fprintf('MN = %.4f, MD = %.4f, SD = %.4f, SE = %.4f, n=%.0f\n',...
      rtdiff_mn{i}, rtdiff_md{i},rtdiff_sd{i},rtdiff_sd{i}./sqrt(length(rtdiff_sd{i})), length(rtdata{i}(:,1)) );

   data_median = rtdiff_md{i};

   % Now calculate what we would expect if the values were randomly assigned.
   nreps = 1000;
   rand_median = zeros(1,nreps);
   rand_mean = zeros(1,nreps);
   total_data = [larger(:)' smaller(:)'];
   for ii = 1:nreps
      index1 = randperm(length(total_data),length(larger) );
      index2 = setdiff(1:length(total_data), index1);
      temp1 = total_data(index1);
      temp2 = total_data(index2);

%       index1 = randperm(length(larger));
%       index2 = randperm(length(larger));
%       temp1 = larger(index1);
%       temp2 = smaller(index2);

      randdiff = abs( temp1 - temp2 );
      rand_median(ii) = nanmedian(randdiff);
      rand_mean(ii) = nanmean(randdiff);
   end % (for i)

   index = find(rand_median < data_median);
   pval = length(index) / nreps;

   subplot(nrows,2,(i)*2);
   hold on;
   [n,x] = hist(rand_median,50);
   edges_rand = center2edge(x);
   n = histc(rand_median, edges_rand);
   n = n ./ sum(n);
   hb = bar(edges_rand,n,'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
   set(hb, 'edgecolor', 'none');
   tickpref;
   minmin = min([data_median min(rand_median)]);
   maxmax = max([data_median max(rand_median)]);
   range = maxmax - minmin;
   xlim([minmin-0.05*range maxmax+0.05*range]);
   ylim([0 max(n)]);
   plot([data_median data_median], ylim, 'k-', 'linewidth', 2);
   xlabel('Median RT Diff');
   title(sprintf('pval = %.4f, n = %.0f', pval, nreps));

end % (for i)

set(gcf,'Position',[114   340   438   365])
print_mfilename(mfilename);

return;
