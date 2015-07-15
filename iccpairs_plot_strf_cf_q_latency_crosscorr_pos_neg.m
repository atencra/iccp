function iccpairs_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs) 
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% iccpairs_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs)
% -----------------------------------------------------------------------
%
%
% Data file: iccpairs-ccpairs.mat
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


cf1 = [ccpairs.cf1];
cf2 = [ccpairs.cf2];

q1 = [ccpairs.q1];
q2 = [ccpairs.q2];

latency1 = [ccpairs.latency1];
latency2 = [ccpairs.latency2];

cf = [cf1(:) cf2(:)];
q = [q1(:) q2(:)];
latency = [latency1(:) latency2(:)];

ccpairs_cf = cf;
ccpairs_q = q;
ccpairs_latency = latency;

if ( min(min(cf)) > 40 )
   cf = cf ./ 1000;
   ccpairs_cf = ccpairs_cf ./ 1000;
end

xmin = [0.4 0 0];
xmax = [20  5    20];
ymin = xmin;
ymax = xmax;
% ymin = [.4  0.25 4];
% ymax = [20  8    17];
xscale = {'log', 'linear', 'linear'};
yscale = {'log', 'linear', 'linear'};

edges{1} = linspace(0, 2, 11);
edges{2} = linspace(0, 4, 11);
edges{3} = linspace(0, 10, 11);

tick = {[0.001 0.01 0.1 1 3 5 10 20], [0:5], [0:4:20]};
ticklabel = {[0.001 0.01 0.1 1 3 5 10 20], [0:5], [0:4:20]};

tick = {[0.625 1.25 2.5 5 10 20], [0:5], [0:4:20]};
ticklabel = {[0.625 1.25 2.5 5 10 20], [0:5], [0:4:20]};

xlabel1 = {'BF (kHz)', 'Q', 'Latency (ms)'};
xlabel2 = {'BF Difference (oct)', 'Q Difference (oct)', 'Latency Difference (ms)'};

spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'cf', 'q', 'latency'};
celltype = {'ccpairs'};

ydiffmax = [0.8 0.5 0.8];
ytick = {linspace(0,ydiffmax(1),5), linspace(0,ydiffmax(2),5), linspace(0,ydiffmax(3),5)}; 
cmap = brewmaps('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 3;
label_left = {'A', 'C', 'E'};
label_right = {'B', 'D', 'F'};

figure;

for i = 1:length(datatype)

   dataspktype = eval(['ccpairs_' datatype{i}]);

   if ( i > 2 )
      dataspktype(:,1) = jitter( dataspktype(:,1) );
      dataspktype(:,2) = jitter( dataspktype(:,2) );
   end 

   dataspktype_pos = dataspktype(index_pos_only,:);
   dataspktype_neg = dataspktype(index_neg_only,:);
   dataspktype_pos_neg = dataspktype(index_pos_neg,:);


   [larger,smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
    [r, p] = corrcoef(log10(larger), log10(smaller));
r
p
length(larger)
pause
   datadiff = abs( larger - smaller );

   [larger_pos,smaller_pos] = vs_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));
   datadiff_pos = abs( larger_pos - smaller_pos );

   [larger_neg,smaller_neg] = vs_largersmaller(dataspktype_neg(:,1),dataspktype_neg(:,2));
   datadiff_neg = abs( larger_neg - smaller_neg );

   [larger_pos_neg,smaller_pos_neg] = ...
        vs_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));
   datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );


   if ( strcmp(datatype{i}, 'cf') )
      datadiff = abs( log2( larger ./ smaller ) );
      datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
      datadiff_neg = abs( log2( larger_neg ./ smaller_neg ) );
      datadiff_pos_neg = abs( log2( larger_pos_neg ./ smaller_pos_neg ) );
%       fprintf('%s  CF r = %.3f, p = %.4f\n', spktype_label{1}, r(1,2), p(1,2) );
   elseif ( strcmp(datatype{i}, 'latency') )
      datadiff = abs( larger - smaller );
      datadiff_pos = abs( larger_pos - smaller_pos );
      datadiff_neg = abs( larger_neg - smaller_neg );
      datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );
%       [r, p] = corrcoef((dataspktype(:,1)), (dataspktype(:,2)) );
%       fprintf('%s  Latency r = %.3f, p = %.4f\n', spktype_label{1}, r(1,2), p(1,2) );
   else % it must be 'q'
      datadiff = abs( log2( larger ./ smaller ) );
      datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
      datadiff_neg = abs( log2( larger_neg ./ smaller_neg ) );
      datadiff_pos_neg = abs( log2( larger_pos_neg ./ smaller_pos_neg ) );
%       [r, p] = corrcoef(log10(dataspktype(:,1)), log10(dataspktype(:,2)) );
%       fprintf('%s  Q r = %.3f, p = %.4f\n', spktype_label{1}, r(1,2), p(1,2) );
   end

   iccpairs_ccc_pos_neg_summary(datadiff_pos, datadiff_neg, datadiff_pos_neg, datatype{i});



   subplot(3,2,(i-1)*2+1);
   hold on;
   plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');

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
   xlabel(sprintf('%s %s', xlabel1{i}, spktype_xlabel{1}));
   ylabel(sprintf('%s %s', xlabel1{i}, spktype_xlabel{1}));
   subplot_label(gca, label_left{i});



   subplot(3,2,2*i);
   hold on;
   n = histc(datadiff_pos, edges{i});
   pdf_pos = n ./ sum(n);
   hp = plot(edges{i}, pdf_pos, 'o-', 'markersize', 3, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
   set(hp, 'color', cmap(1,:));

   n = histc(datadiff_neg, edges{i});
   pdf_neg = n ./ sum(n);
   hp = plot(edges{i}, pdf_neg, 'o-', 'markersize', 3, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
   set(hp, 'color', cmap(2,:));

   n = histc(datadiff_pos_neg, edges{i});
   pdf_pos_neg = n ./ sum(n);
   hp = plot(edges{i}, pdf_pos_neg, 'o-', 'markersize', 3, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
   set(hp, 'color', cmap(3,:));

   box off;
   tickpref;
   xtick = edges{i}(1:2:end);
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   range = max(edges{i}) - min(edges{i});
   xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
   ylim([0 ydiffmax(i)]);
   ylabel('Proportion');
   xlabel(xlabel2{i});
   if ( i == 1 ), legend('Pos', 'Neg', 'Pos+Neg'); end;
   subplot_label(gca,label_right{i});


end % for i

set(gcf,'position', [1132 402 438 517]);
print_mfilename(mfilename);


return;






% 
% 
% 
% 
% function plot_pairs_population_cf_q_latency(cf, bw, q, latency) 
% % plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
% %
% % plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% % -----------------------------------------------------------------------
% %
% %
% %
% % caa 4/21/11
% 
% index = find(cf(:,1)>500 & cf(:,2)>500 );
% cf = cf(index,:);
% q = q(index,:);
% bw = bw(index,:);
% latency = latency(index,:);
% 
% 
% % cf
% datastr(1).data = cf;
% datastr(1).datatype = 'cf';
% datastr(1).xmin = 500;
% datastr(1).xmax = 26000;
% datastr(1).ymin = 500;
% datastr(1).ymax = 26000;
% datastr(1).ydiffmax = 0.75;
% datastr(1).ydifftick = 0:0.25:0.75;
% datastr(1).xscale = 'log';
% datastr(1).yscale = 'log';
% datastr(1).edges = linspace(0, 2.8, 15);
% datastr(1).tick = 500*2.^(0:1:5);
% datastr(1).ticklabel = 500*2.^(0:1:5)./1000;
% datastr(1).xlabel1 = 'BF (Hz)';
% datastr(1).xlabel2 = 'BF Diff (oct)';
% 
% % bw
% datastr(2).data = bw;
% datastr(2).datatype = 'bw';
% datastr(2).xmin = 0;
% datastr(2).xmax = 8;
% datastr(2).ymin = 0;
% datastr(2).ymax = 8;
% datastr(2).ydiffmax = 0.305;
% datastr(2).ydifftick = 0:0.1:0.3;
% datastr(2).xscale = 'linear';
% datastr(2).yscale = 'linear';
% datastr(2).edges = linspace(0, 4, 17);
% datastr(2).tick = [0:8];
% datastr(2).ticklabel = [0:8];
% datastr(2).xlabel1 = 'BW (oct)';
% datastr(2).xlabel2 = 'BW Diff (oct)';
% 
% 
% % q
% datastr(3).data = q;
% datastr(3).datatype = 'q';
% datastr(3).xmin = 0;
% datastr(3).xmax = 5;
% datastr(3).ymin = 0;
% datastr(3).ymax = 5;
% datastr(3).ydiffmax = 0.3;
% datastr(3).ydifftick = 0:0.1:0.3;
% datastr(3).xscale = 'linear';
% datastr(3).yscale = 'linear';
% datastr(3).edges = linspace(0, 4, 17);
% datastr(3).tick = [0:5];
% datastr(3).ticklabel = [0:5];
% datastr(3).xlabel1 = 'Q';
% datastr(3).xlabel2 = 'Q Diff (oct)';
% 
% 
% % latency
% datastr(4).data = latency;
% datastr(4).datatype = 'latency';
% datastr(4).xmin = 0;
% datastr(4).xmax = 25;
% datastr(4).ymin = 0;
% datastr(4).ymax = 25;
% datastr(4).ydiffmax = 0.8;
% datastr(4).ydifftick = 0:0.2:0.8;
% datastr(4).xscale = 'linear';
% datastr(4).yscale = 'linear';
% datastr(4).edges = linspace(0, 20, 17);
% datastr(4).tick = [0:5:25];
% datastr(4).ticklabel = [0:5:25];
% datastr(4).xlabel1 = 'Latency (ms)';
% datastr(4).xlabel2 = 'Latency Diff (ms)';
% 
% 
% figure;
% nrows = length(datastr);
% for i = 1:length(datastr)
% 
%    data = datastr(i).data;
%    datatype = datastr(i).datatype;
%    xmin = datastr(i).xmin;
%    xmax = datastr(i).xmax;
%    ymin = datastr(i).ymin;
%    ymax = datastr(i).ymax;
%    ydiffmax = datastr(i).ydiffmax;
%    ydifftick = datastr(i).ydifftick;
%    xscale = datastr(i).xscale;
%    yscale = datastr(i).yscale;
%    edges = datastr(i).edges;
%    tick = datastr(i).tick;
%    ticklabel = datastr(i).ticklabel;
%    xlabel1 = datastr(i).xlabel1;
%    xlabel2 = datastr(i).xlabel2;
% 
%    if ( strcmp(datatype, 'cf') )
%       datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
%    elseif ( strcmp(datatype, 'latency') )
%       datadiff = abs( data(:,1) - data(:,2) );
%    else % it must be 'q'
%       datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
% %       datadiff = abs( data(:,1) - data(:,2) );
%    end
% 
% 
%    subplot(nrows,2,(i-1)*2+1);
%    hold on;
%    [larger, smaller] = vs_largersmaller(data(:,1), data(:,2));
%    plot(larger, smaller, 'ko', 'markerfacecolor', 'k', 'markersize', 1);
%    plot([xmin xmax], [ymin ymax], 'k-');
%    set(gca,'xscale', xscale, 'yscale', yscale);
%    tickpref;
%    set(gca,'xtick', tick, 'xticklabel', ticklabel);
%    set(gca,'ytick', tick, 'yticklabel', ticklabel);
%    xlim([xmin xmax]);
%    ylim([ymin ymax]);
%    xlabel(sprintf('%s Neuron 1', xlabel1));
%    ylabel(sprintf('%s Neuron 2', xlabel1));
% 
% 
%    subplot(nrows,2,2*i);
%    n = histc(datadiff, edges);
%    n = n ./ sum(n);
%    hb = bar(edges, n, 'histc');
%    set(hb, 'facecolor', [0.7 0.7 0.7]);
%    box off;
%    tickpref;
%    xtick = edges(1:2:end);
%    set(gca,'xtick', xtick, 'xticklabel', xtick);
%    ytick = ydifftick;
%    set(gca,'ytick', ytick, 'yticklabel', ytick);
%    range = max(edges) - min(edges);
%    xlim([min(edges) max(edges)]);
%    ylim([0 ydiffmax]);
%    ylabel('Proportion');
%    xlabel(xlabel2);
% 
%    fprintf('\n');
%    fprintf('%s\n', xlabel1);
%    [r,p] = corrcoef(larger, smaller);
%    fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
%    fprintf('%s\n', xlabel2);
%    fprintf('N=%.0f,MN=%.4f,SD=%.4f', size(data,1),mean(datadiff), std(datadiff));
%    fprintf('MD=%.4f,MAD=%.4f,MX=%.4f', median(datadiff), mad(datadiff), max(datadiff));
%    fprintf('\n\n');
% 
% end % for i
% 
% set(gcf,'position',[996   229   438   749]);
% print_mfilename(mfilename);
% 
% return;
% 
% 
