function iccpairs_plot_ptdata_cf_q_latency(ptdata) 
% iccpairs_plot_ptdata_bw ICC Pairwise FRA Q comparison
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




cf = ptdata.cf;
q = ptdata.q;
bw = ptdata.bw;
latency = ptdata.latency;

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
datastr(1).ydiffmax = 0.6;
datastr(1).ydifftick = 0:0.2:0.6;
datastr(1).xscale = 'log';
datastr(1).yscale = 'log';
datastr(1).edges = linspace(0, 1.2, 13);
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
datastr(2).ydiffmax = 0.4;
datastr(2).ydifftick = 0:0.1:0.4;
datastr(2).xscale = 'linear';
datastr(2).yscale = 'linear';
datastr(2).edges = linspace(0, 3.5, 15);
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
datastr(3).ydiffmax = 0.4;
datastr(3).ydifftick = 0:0.1:0.4;

datastr(3).xscale = 'linear';
datastr(3).yscale = 'linear';
datastr(3).edges = linspace(0, 3.5, 15);
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
datastr(4).edges = linspace(0, 15, 16);
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
   plot(larger, smaller, 'ko', 'markerfacecolor', [0.6 0.6 0.6], 'markersize', 2);
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
   xtick = edges(1:3:end);
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

set(gcf,'position',[996   229   361 551]);
print_mfilename(mfilename);

return;

