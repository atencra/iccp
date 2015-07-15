function plot_allpairs_spktype_crosscorr_params(ccpairs)
%plot_pairs_spktype_crosscorr_params  Values from same channel cross-corr functions
%
% plot_pairs_spktype_crosscorr_params(ccpairs)
%
% Modified from plot_pairs_spktype_crosscorr_params to work just with
% ccpairs data, and not fsufsu, fsursu, or rsursu data.

% Functions called:
% cschemes.m
% splot.m
% smooth3p.m
% edge2center.m
% rotateticklabel.m

% Input:
% ----------------------
% 1) ccpairs: ccpairs data loaded from icc_crosscorr_strf_new2

% Output:
% ----------------------
% 1) Plots hisograms of Peak Delay, Half Width, and Correlation Coefficient
% 2) Plots mean values for the above parameters as well.

% Preliminary Values
ccpairs_pos = [ccpairs.position];
% ccpairs_rho = [ccpairs.rho];
ccpairs_rho = [ccpairs.ccc];
ccpairs_rho(ccpairs_rho<=0) = 0.0001;
ccpairs_pd = [ccpairs.peakdelay];
ccpairs_pd = abs(ccpairs_pd);
ccpairs_hw = [ccpairs.halfwidth];
ccpairs_sig = logical([ccpairs.significant]);

% [min(ccpairs_pd) max(ccpairs_pd)]
% length(ccpairs_pd)
% figure;
% hist(ccpairs_pd)
% pause

ccpairs_pos_sig = ccpairs_pos(ccpairs_sig);
ccpairs_rho_sig = ccpairs_rho(ccpairs_sig);
ccpairs_pd_sig = ccpairs_pd(ccpairs_sig);
ccpairs_hw_sig = ccpairs_hw(ccpairs_sig);


% Find monosynaptic FSU-FSU pairs
monoindex = find( abs(ccpairs_pd_sig)>0 & abs(ccpairs_pd_sig<5) & ...
   (ccpairs_hw_sig<=10) );
ccpairs_rho_mono = ccpairs_rho_sig(monoindex);
ccpairs_pd_mono = ccpairs_pd_sig(monoindex);
ccpairs_hw_mono = ccpairs_hw_sig(monoindex);
ccpairs_pos_mono = ccpairs_pos_sig(monoindex);


fprintf('\n');
fprintf('CCPairs   Mono  #Sig  Total\n');
fprintf('          %4.0f  %4.0f  %4.0f\n', ...
length(ccpairs_pos_mono), length(ccpairs_pos_sig), length(ccpairs_pos));
fprintf('\n');

close all;

splotnum = [4 2; 10 8; 16 14];
edgescell = {linspace(0,5,11), linspace(0,8,17), ...
   10.^( linspace(log10(0.001), log10(0.6), 19)) };
xtick = {linspace(0,5,6), linspace(0,8,5)};
ytick = {linspace(0,0.6,4), linspace(0,0.6,4),linspace(0,0.2,5)};
xlabelstr = {'Peak Delay (ms)', 'Half Width (ms)', 'Corr Coef'};

cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 's', '^'};
markersize = 4;
linewidth = 2;
ymax = [0.65 0.65 0.2];
datatype = {'pd', 'hw', 'rho'};

for i = 1:3

%    fsufsu_data = abs(eval(['fsufsu_' datatype{i} '_sig']));
   ccpairs_data = abs(eval(['ccpairs_' datatype{i}]));

%    celltypestr = {'fsufsu_data', 'fsursu_data', 'rsursu_data'};
   celltypestr = {'ccpairs_data'};

   splot(6,3,splotnum(i,:));
   edges = edgescell{i};
   hold on;

   for j = 1:length(celltypestr)
      celltypedata = eval(celltypestr{j});

      n = histc(celltypedata, edges);
      n = n ./ sum(n);
      if j == length(celltypestr)
         n = smooth3p(n);
      end
      centers = edge2center(edges);
      plot(centers, n(1:end-1), 'o-', 'color', cmap(j,:), ...
         'markersize', markersize, 'markerfacecolor', cmap(j,:), ...
         'markeredgecolor', cmap(j,:), 'linewidth', linewidth);
   end % (for j)

   box off;
   tickpref;
   if ( i == 1 || i == 2 )
      set(gca,'xtick', xtick{i}, 'xticklabel', xtick{i});
   end
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});

   if ( i == 1 || i == 2 )
      range = max(edges) - min(edges);
      xlim([min(edges)-0.025*range max(edges)+0.025*range]);
   end

   if ( i == 3 )
      xlim([min(edges) max(edges)]);
      set(gca,'xscale', 'log');
   end

   ylim([0 ymax(i)]);
   xlabel(xlabelstr{i});
   ylabel('Proportion');

   if ( i == 1 )
%       legend('FSU-FSU', 'FSU-RSU', 'RSU-RSU');
%       legend('FSU-RSU', 'RSU-RSU');
   end

end % (for i )




datalabel = {'Peak Delay (ms)', 'Half Width (ms)', 'Correlation Coefficient'};
datatype = {'pd', 'hw', 'rho'};
ydatalabel = {'Mean PD (ms)', 'Mean HW (ms)', 'Mean CC'};

ymax = [4.2 3 0.07];
ytickmax = [4 3 0.06];
plotnum = [3 9 15];

for i = 1:3

%    fsufsu_data = abs(eval(['fsufsu_' datatype{i} '_sig']));
   ccpairs_data = abs(eval(['ccpairs_' datatype{i}]));

   ccpairs_mn = nanmean(ccpairs_data);
   

   ccpairs_sd = nanstd(ccpairs_data);
 

   ccpairs_se = ccpairs_sd / sqrt(length(ccpairs_data));
   
   fprintf('\n%s:\n', datalabel{i});
   fprintf('CCPairs mean (SD): %.3f (%.3f)\n', ccpairs_mn, ccpairs_sd);
   p = ranksum(log10(ccpairs_data(~isnan(ccpairs_data))), log10(ccpairs_data));
   fprintf('Neuron 1 vs Neuron 2: p = %.4f\n', p);


   splot(6,3,plotnum(i));
   hold on;
   hb = bar([0:2], [0 ccpairs_mn 0]);
   set(hb, 'facecolor', [0.7 0.7 0.7]);

   plot([1 1], [ccpairs_mn ccpairs_mn+ccpairs_se], 'k-');
   plot([0.8 1.2], [ccpairs_mn+ccpairs_se ccpairs_mn+ccpairs_se], 'k-');

% %     What do I do for these
%    plot([2 2], [fsursu_mn fsursu_mn+fsursu_se], 'k-');
%    plot([1.8 2.2], [fsursu_mn+fsursu_se fsursu_mn+fsursu_se], 'k-');
% 
%    plot([3 3], [rsursu_mn rsursu_mn+rsursu_se], 'k-');
%    plot([2.8 3.2], [rsursu_mn+rsursu_se rsursu_mn+rsursu_se], 'k-');

%    hb = bar([0:3], [0 fsursu_mn rsursu_mn 0]);
%    set(hb, 'facecolor', [0.7 0.7 0.7]);
%    plot([2 2]-1, [fsursu_mn fsursu_mn+fsursu_se], 'k-');
%    plot([1.8 2.2]-1, [fsursu_mn+fsursu_se fsursu_mn+fsursu_se], 'k-');
%    plot([3 3]-1, [rsursu_mn rsursu_mn+rsursu_se], 'k-');
%    plot([2.8 3.2]-1, [rsursu_mn+rsursu_se rsursu_mn+rsursu_se], 'k-');

   box off;
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   xtick = 1;
   xticklabel = {'CCPairs'};
%    xticklabel = {'F-R', 'R-R'};
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   rotateticklabel(gca,90);
   ytick = linspace(0,ytickmax(i),3);
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   xlim([0 4]);
   ylim([0 ymax(i)]);
   ylabel(ydatalabel{i});
   title(datalabel{i});

end % (for i)


set(gcf,'position', [55   330   440   660]);
print_mfilename(mfilename);






close all;

edgescell = {0:1.5:20, linspace(0,8,17), ...
   10.^( linspace(log10(0.001), log10(0.6), 19)) };
xtick = {0:3:20, linspace(0,8,5)};
ytick = {linspace(0,0.4,5), linspace(0,0.4,5),linspace(0,0.2,5)};
xlabelstr = {'Peak Delay (ms)', 'Half Width (ms)', 'Corr Coef'};

cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 's', '^'};
markersize = 4;
linewidth = 2;
ymax = [0.4 0.45 0.2];
datatype = {'pd', 'hw', 'rho'};

for i = 1:3

   ccpairs_data = abs(eval(['ccpairs_' datatype{i}]));
   ccpairs_mn = nanmean(ccpairs_data);
   ccpairs_md = nanmedian(ccpairs_data);
   ccpairs_mad = mad(ccpairs_data(~isnan(ccpairs_data)));
   ccpairs_sd = nanstd(ccpairs_data);
   ccpairs_se = ccpairs_sd / sqrt(length(ccpairs_data));

   edges = edgescell{i};
   n = histc(ccpairs_data, edges);
   n = n ./ sum(n);
%          n = smooth3p(n);
   centers = edge2center(edges);

   subplot(3,1,i);
   hold on;
   hb = bar(edges,n, 'histc');
   set(hb, 'facecolor', 0.6 * ones(1,3));
   plot([ccpairs_mn ccpairs_mn], [0 ymax(i)], 'k');
   box off;
   tickpref;
   if ( i == 1 || i == 2 )
      set(gca,'xtick', xtick{i}, 'xticklabel', xtick{i});
   end
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   xlim([min(edges) max(edges)]);

   if ( i == 1 || i == 2 )
      xlim([min(edges) max(edges)]);
   end

   if ( i == 3 )
      xlim([min(edges) max(edges)]);
      set(gca,'xscale', 'log');
%       set(hb, 'markersize', 0)
      h = findobj(gca,'Type','line');
      set(h,'Marker','none'); 
   end

   if ( i == 1 )
      text(6, 0.3, sprintf('mn=%.2f, sd=%.2f\nse=%.2f md=%.2f, mad=%.2f',...
         ccpairs_mn, ccpairs_sd, ccpairs_se, ccpairs_md, ccpairs_mad));
   end

   if ( i == 2 )
      text(2, 0.4, sprintf('mn=%.2f, sd=%.2f\nse=%.2f md=%.2f, mad=%.2f',...
         ccpairs_mn, ccpairs_sd, ccpairs_se, ccpairs_md, ccpairs_mad));
   end

   if ( i == 3 )
      text(0.002, 0.17, sprintf('mn=%.2f, sd=%.2f\nse=%.2f md=%.2f, mad=%.2f',...
         ccpairs_mn, ccpairs_sd, ccpairs_se, ccpairs_md, ccpairs_mad));
   end

   ylim([0 ymax(i)]);
   xlabel(xlabelstr{i});
   ylabel('Proportion');

end % (for i )

set(gcf,'position', [1100 300 321 660]);
print_mfilename(mfilename);




return;




