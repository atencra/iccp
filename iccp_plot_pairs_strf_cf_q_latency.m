function [cf, bw, q, latency] = iccp_plot_pairs_strf_cf_q_latency(cf, bw, q, latency)
% iccp_plot_pairs_strf_cf_q_latency - compare pure tone params for pairs of neurons
%
% 
% iccp_plot_pairs_strf_cf_q_latency
% --------------------------------------------------------
% The function finds the neurons that were recorded from the 
% same channel, and then compares the pure tone parameters
% for these neurons.
%
% The analysis tells us the variability of the parameters in
% a locally confined area of cortex.
%
% plot_pairs_strf_cf_q_latency;



% BF (Hz)
% r=0.926
% BF Diff (oct)
% N=1045,MN=0.20,SD=0.29MD=0.10,MAD=0.08,MX=2.61
% 
% BW (oct)
% r=0.542
% BW Diff (oct)
% N=1045,MN=0.68,SD=0.60MD=0.53,MAD=0.37,MX=3.88
% 
% Q
% r=0.691
% Q Diff (oct)
% N=1045,MN=0.73,SD=0.64MD=0.58,MAD=0.38,MX=4.44
% 
% Latency (ms)
% r=0.845
% Latency Diff (ms)
% N=1045,MN=1.13,SD=1.68MD=0.50,MAD=0.50,MX=16.00








close all;

if ( nargin ~= 4 )
   [position, cf, bw, q, latency] = get_pairs_spktype_strf_cf_q_latency;
end



% [position, cf, q, latency, dur1, dur2] = get_pairs_spktype_strf_cf_q_latency;

% [cfdata, qdata, latdata] = plot_pairs_population_diff_layers_cf_q_latency(position, cf, q, latency);

% [cfdata, qdata, latdata] = plot_pairs_cdf_diff_pop_layers_cf_q_latency(position, cf, q, latency);



% close all;

% plot_pairs_spktype_population_cf_q_latency(cf, q, latency, dur1, dur2);

% plot_pairs_spktype_fsu_rsu_cdf_diff_cf_q_latency(cf, q, latency, dur1, dur2);



plot_pairs_population_cf_q_latency(cf, bw, q, latency); 

plot_pairs_population_cf_q_latency_diff_rand(cf, bw, q, latency);


% plot_pairs_layer_cf_q_latency(position, cf, 'cf');
% plot_pairs_layer_cf_q_latency(position, latency, 'latency');
% plot_pairs_layer_cf_q_latency(position, q, 'q');

return;



function [position, cf, bw, q, latency] = get_pairs_spktype_strf_cf_q_latency

cf = [];
bw = [];
q = [];
latency = [];
position = [];
% dur1 = [];
% dur2 = [];

d = dir('*-strfcmb-pairs-ptparams.mat');

for n = 1:length(d)

	filename = d(n).name;
	s = load(filename, 'ptparams');
   ptparams = s.ptparams;

%    index = findstr(filename, '-ptparams');
%    prefix = filename(1:index-1);

%    dspktype = dir(sprintf('%s-spktype.mat',prefix));
%    filename = dspktype.name;
%    s = load(filename, 'spktype');
%    spktype = s.spktype;
% 
%    if ( length(spktype) ~= length(ptparams) )
%       error('spktype and ptparams are not the same length.');
%    end


   chan = [ptparams.chan];
   chan_unique = unique(chan);

   for i = 1:length(chan_unique)

      index = find( chan_unique(i) == chan );

      if ( length(index) > 1 )

         cmb = nchoosek(index, 2); % determine all possible pairs 

         [nr, nc] = size(cmb);

         for j = 1:nr

            index1 = cmb(j,1);
            index2 = cmb(j,2);

            % First get the spike shape data

            % Neuron #1, phase 1 and 2
   %          phase11 = spktype(index1).zero_crossings(2) - spktype(index1).zero_crossings(1);
   %          phase12 = spktype(index1).zero_crossings(3) - spktype(index1).zero_crossings(2);

            % Neuron #2, phase 1 and 2
   %          phase21 = spktype(index2).zero_crossings(2) - spktype(index2).zero_crossings(1);
   %          phase22 = spktype(index2).zero_crossings(3) - spktype(index2).zero_crossings(2);


   %          dur1 = [dur1; phase11 + phase12];
   %          dur2 = [dur2; phase21 + phase22];

   %          chan1 = spktype(index1).chan(1); 
   %          chan2 = spktype(index2).chan(1); 
   %          model1 = spktype(index1).model(1); 
   %          model2 = spktype(index2).model(1); 

            chan3 = ptparams(index1).chan(1); 
            chan4 = ptparams(index2).chan(1); 
            model3 = ptparams(index1).model(1); 
            model4 = ptparams(index2).model(1); 

            pos = ptparams(index1).position;
            position = [position; pos];

            cf = [cf; ptparams(index1).cf ptparams(index2).cf];
            q = [q; ptparams(index1).q ptparams(index2).q];

            bw1 = abs( log2(ptparams(index1).fupper ./ ptparams(index1).flower  ) );
            bw2 = abs( log2(ptparams(index2).fupper ./ ptparams(index2).flower  ) );
            bw = [bw; bw1 bw2];

            latency1 = ptparams(index1).latency;
            latency2 = ptparams(index2).latency;

            latency = [latency; latency1 latency2];

         end % (for j)

      end % (if)

   end % (for i)

end % (for n)

% Now sort the data so that dur1 has the lowest durations, dur2 the
% longest. Swap the data values accordingly

% for i = 1:length(dur1)
% 
%    if ( dur1(i) > dur2(i) )
% 
%       a = dur1(i);
%       b = dur2(i);
%       dur1(i) = b;
%       dur2(i) = a;
% 
%       a = cf(i,1);
%       b = cf(i,2);
%       cf(i,1) = b;
%       cf(i,2) = a;
%    
%       if ( max(cf(:)) > 1000 )
%          cf = cf ./ 1000;
%       end
% 
%       a = q(i,1);
%       b = q(i,2);
%       q(i,1) = b;
%       q(i,2) = a;
% 
%       a = latency(i,1);
%       b = latency(i,2);
%       latency(i,1) = b;
%       latency(i,2) = a;
% 
%    end
% 
% end % (for i)

return;




function plot_pairs_population_cf_q_latency(cf, bw, q, latency) 
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% -----------------------------------------------------------------------
%
%
%
% caa 4/21/11

index = find(cf(:,1)>500 & cf(:,2)>500 );
cf = cf(index,:);
q = q(index,:);
bw = bw(index,:);
latency = latency(index,:);


% cf
datastr(1).data = cf;
datastr(1).datatype = 'cf';
datastr(1).xmin = 500;
datastr(1).xmax = 26000;
datastr(1).ymin = 500;
datastr(1).ymax = 26000;
datastr(1).ydiffmax = 0.75;
datastr(1).ydifftick = 0:0.25:0.75;
datastr(1).xscale = 'log';
datastr(1).yscale = 'log';
datastr(1).edges = linspace(0, 2.8, 15);
datastr(1).tick = 500*2.^(0:1:5);
datastr(1).ticklabel = 500*2.^(0:1:5)./1000;
datastr(1).xlabel1 = 'BF (Hz)';
datastr(1).xlabel2 = 'BF Diff (oct)';

% bw
datastr(2).data = bw;
datastr(2).datatype = 'bw';
datastr(2).xmin = 0;
datastr(2).xmax = 8;
datastr(2).ymin = 0;
datastr(2).ymax = 8;
datastr(2).ydiffmax = 0.305;
datastr(2).ydifftick = 0:0.1:0.3;
datastr(2).xscale = 'linear';
datastr(2).yscale = 'linear';
datastr(2).edges = linspace(0, 4, 17);
datastr(2).tick = [0:8];
datastr(2).ticklabel = [0:8];
datastr(2).xlabel1 = 'BW (oct)';
datastr(2).xlabel2 = 'BW Diff (oct)';


% q
datastr(3).data = q;
datastr(3).datatype = 'q';
datastr(3).xmin = 0;
datastr(3).xmax = 5;
datastr(3).ymin = 0;
datastr(3).ymax = 5;
datastr(3).ydiffmax = 0.3;
datastr(3).ydifftick = 0:0.1:0.3;
datastr(3).xscale = 'linear';
datastr(3).yscale = 'linear';
datastr(3).edges = linspace(0, 4, 17);
datastr(3).tick = [0:5];
datastr(3).ticklabel = [0:5];
datastr(3).xlabel1 = 'Q';
datastr(3).xlabel2 = 'Q Diff (oct)';


% latency
datastr(4).data = latency;
datastr(4).datatype = 'latency';
datastr(4).xmin = 0;
datastr(4).xmax = 25;
datastr(4).ymin = 0;
datastr(4).ymax = 25;
datastr(4).ydiffmax = 0.8;
datastr(4).ydifftick = 0:0.2:0.8;
datastr(4).xscale = 'linear';
datastr(4).yscale = 'linear';
datastr(4).edges = linspace(0, 20, 17);
datastr(4).tick = [0:5:25];
datastr(4).ticklabel = [0:5:25];
datastr(4).xlabel1 = 'Latency (ms)';
datastr(4).xlabel2 = 'Latency Diff (ms)';


figure;
nrows = length(datastr);
for i = 1:length(datastr)

   data = datastr(i).data;
   datatype = datastr(i).datatype;
   xmin = datastr(i).xmin;
   xmax = datastr(i).xmax;
   ymin = datastr(i).ymin;
   ymax = datastr(i).ymax;
   ydiffmax = datastr(i).ydiffmax;
   ydifftick = datastr(i).ydifftick;
   xscale = datastr(i).xscale;
   yscale = datastr(i).yscale;
   edges = datastr(i).edges;
   tick = datastr(i).tick;
   ticklabel = datastr(i).ticklabel;
   xlabel1 = datastr(i).xlabel1;
   xlabel2 = datastr(i).xlabel2;

   if ( strcmp(datatype, 'cf') )
      datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
   elseif ( strcmp(datatype, 'latency') )
      datadiff = abs( data(:,1) - data(:,2) );
   else % it must be 'q'
      datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
%       datadiff = abs( data(:,1) - data(:,2) );
   end


   subplot(nrows,2,(i-1)*2+1);
   hold on;
   [larger, smaller] = vs_largersmaller(data(:,1), data(:,2));
   plot(larger, smaller, 'ko', 'markerfacecolor', 'k', 'markersize', 1);
   plot([xmin xmax], [ymin ymax], 'k-');
   set(gca,'xscale', xscale, 'yscale', yscale);
   tickpref;
   set(gca,'xtick', tick, 'xticklabel', ticklabel);
   set(gca,'ytick', tick, 'yticklabel', ticklabel);
   xlim([xmin xmax]);
   ylim([ymin ymax]);
   xlabel(sprintf('%s Neuron 1', xlabel1));
   ylabel(sprintf('%s Neuron 2', xlabel1));


   subplot(nrows,2,2*i);
   n = histc(datadiff, edges);
   n = n ./ sum(n);
   hb = bar(edges, n, 'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
   box off;
   tickpref;
   xtick = edges(1:2:end);
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   ytick = ydifftick;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   range = max(edges) - min(edges);
   xlim([min(edges) max(edges)]);
   ylim([0 ydiffmax]);
   ylabel('Proportion');
   xlabel(xlabel2);

   fprintf('\n');
   fprintf('%s\n', xlabel1);
   [r,p] = corrcoef(larger, smaller);
   fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
   fprintf('%s\n', xlabel2);
   fprintf('N=%.0f,MN=%.4f,SD=%.4f', size(data,1),mean(datadiff), std(datadiff));
   fprintf('MD=%.4f,MAD=%.4f,MX=%.4f', median(datadiff), mad(datadiff), max(datadiff));
   fprintf('\n\n');

end % for i

set(gcf,'position',[996   229   438   749]);
print_mfilename(mfilename);

return;




function plot_pairs_population_cf_q_latency_diff_rand(cf, bw, q, latency) 
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% -----------------------------------------------------------------------
%
%
%
% caa 4/21/11

index = find(cf(:,1)>500 & cf(:,2)>500 );
cf = cf(index,:);
q = q(index,:);
bw = bw(index,:);
latency = latency(index,:);


% cf
datastr(1).data = cf;
datastr(1).datatype = 'cf';
datastr(1).xmin = 500;
datastr(1).xmax = 26000;
datastr(1).ymin = 500;
datastr(1).ymax = 26000;
datastr(1).ydiffmax = 0.75;
datastr(1).ydifftick = 0:0.25:0.75;
datastr(1).xscale = 'log';
datastr(1).yscale = 'log';
datastr(1).edges = linspace(0, 2.8, 15);
datastr(1).tick = 500*2.^(0:1:5);
datastr(1).ticklabel = 500*2.^(0:1:5)./1000;
datastr(1).xlabel1 = 'BF (Hz)';
datastr(1).xlabel2 = 'BF Diff (oct)';

% bw
datastr(2).data = bw;
datastr(2).datatype = 'bw';
datastr(2).xmin = 0;
datastr(2).xmax = 8;
datastr(2).ymin = 0;
datastr(2).ymax = 8;
datastr(2).ydiffmax = 0.305;
datastr(2).ydifftick = 0:0.1:0.3;
datastr(2).xscale = 'linear';
datastr(2).yscale = 'linear';
datastr(2).edges = linspace(0, 4, 17);
datastr(2).tick = [0:8];
datastr(2).ticklabel = [0:8];
datastr(2).xlabel1 = 'BW (oct)';
datastr(2).xlabel2 = 'BW Diff (oct)';


% q
datastr(3).data = q;
datastr(3).datatype = 'q';
datastr(3).xmin = 0;
datastr(3).xmax = 5;
datastr(3).ymin = 0;
datastr(3).ymax = 5;
datastr(3).ydiffmax = 0.3;
datastr(3).ydifftick = 0:0.1:0.3;
datastr(3).xscale = 'linear';
datastr(3).yscale = 'linear';
datastr(3).edges = linspace(0, 4, 17);
datastr(3).tick = [0:5];
datastr(3).ticklabel = [0:5];
datastr(3).xlabel1 = 'Q';
datastr(3).xlabel2 = 'Q Diff (oct)';


% latency
datastr(4).data = latency;
datastr(4).datatype = 'latency';
datastr(4).xmin = 0;
datastr(4).xmax = 25;
datastr(4).ymin = 0;
datastr(4).ymax = 25;
datastr(4).ydiffmax = 0.8;
datastr(4).ydifftick = 0:0.2:0.8;
datastr(4).xscale = 'linear';
datastr(4).yscale = 'linear';
datastr(4).edges = linspace(0, 20, 17);
datastr(4).tick = [0:5:25];
datastr(4).ticklabel = [0:5:25];
datastr(4).xlabel1 = 'Latency (ms)';
datastr(4).xlabel2 = 'Latency Diff (ms)';


figure;
nrows = length(datastr);
for i = 1:length(datastr)

   data = datastr(i).data;
   datatype = datastr(i).datatype;
   xmin = datastr(i).xmin;
   xmax = datastr(i).xmax;
   ymin = datastr(i).ymin;
   ymax = datastr(i).ymax;
   ydiffmax = datastr(i).ydiffmax;
   ydifftick = datastr(i).ydifftick;
   xscale = datastr(i).xscale;
   yscale = datastr(i).yscale;
   edges = datastr(i).edges;
   tick = datastr(i).tick;
   ticklabel = datastr(i).ticklabel;
   xlabel1 = datastr(i).xlabel1;
   xlabel2 = datastr(i).xlabel2;

   if ( strcmp(datatype, 'cf') )
      datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
   elseif ( strcmp(datatype, 'latency') )
      datadiff = abs( data(:,1) - data(:,2) );
   else % it must be 'q' or 'bw'
      datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
%       datadiff = abs( data(:,1) - data(:,2) );
   end
   data_median = median(datadiff);

   [larger, smaller] = vs_largersmaller(data(:,1), data(:,2));


   subplot(nrows,2,(i-1)*2+1);
   n = histc(datadiff, edges);
   n = n ./ sum(n);
   hb = bar(edges, n, 'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
   box off;
   tickpref;
   xtick = edges(1:2:end);
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   ytick = ydifftick;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   range = max(edges) - min(edges);
   xlim([min(edges) max(edges)]);
   ylim([0 ydiffmax]);
   ylabel('Proportion');
   xlabel(xlabel2);

   fprintf('\n');
   fprintf('%s\n', xlabel1);
   [r,p] = corrcoef(larger, smaller);
   fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
   fprintf('%s\n', xlabel2);
   fprintf('N=%.0f,MN=%.4f,SD=%.4f', size(data,1),mean(datadiff), std(datadiff));
   fprintf('MD=%.4f,MAD=%.4f,MX=%.4f', median(datadiff), mad(datadiff), max(datadiff));
   fprintf('\n\n');


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
      if ( strcmp(datatype, 'cf') )
         randdiff = abs( log2( temp1 ./ temp2 ) );
      elseif ( strcmp(datatype, 'latency') )
         randdiff = abs( temp1 - temp2 );
      else % it must be 'q' or 'bw'
         randdiff = abs( log2( temp1 ./ temp2 ) );
      end
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
   xlabel(sprintf('Median %s', xlabel2));
   title(sprintf('pval = %.4f, n = %.0f', pval, nreps));

end % for i

set(gcf,'position',[1405 147 438 749]);
print_mfilename(mfilename);

return;







