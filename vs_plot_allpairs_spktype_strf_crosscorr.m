function vs_plot_allpairs_spktype_strf_crosscorr(ccpairs)
% plot_allpairs_spktype_strf_crosscorr STRF params for cross-corr cell types
%
% 
% plot_allpairs_spktype_crosscorr(ccpairs)
% --------------------------------------------------------
%

% Modified from plot_pairs_spkytpe_strf_crosscorr() to work with just
% ccpairs data. Performs the same function as
% plot_pairs_spktype_strf_crosscorr, just with a different input.

% Uses ccpairs data from icc_crosscorr_strf_new2.mat

close all;

% load icc_crosscorr_strf_new3.mat

% iccpairs_plot_strf_similarity_crosscorr_pos_neg(ccpairs);

% iccpairs_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs);

% plot_allpairs_spktype_strf_crosscorr_similarity(ccpairs);

% iccpairs_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs);

% plot_allpairs_spktype_strf_crosscorr_strf_fr_pli(ccpairs);

% plot_allpairs_spktype_strf_crosscorr_mtf_btmf(ccpairs);

% plot_allpairs_spktype_strf_crosscorr_fio_similarity(ccpairs);

% plot_allpairs_spktype_strf_crosscorr_fio_asi(ccpairs);

plot_allpairs_spktype_strf_crosscorr_fio_asi_similarity(ccpairs);

% plot_allpairs_spktype_strf_crosscorr_mtf_btmf_bsmf(ccpairs);


return;



%------------------- Function Definitions ------------------------




function iccpairs_plot_strf_similarity_crosscorr_pos_neg(ccpairs)

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

sum(sigfeature_pos & sigfeature_neg)
sum(sigfeature_pos & ~sigfeature_neg)
sum(~sigfeature_pos & sigfeature_neg)
sum(sigfeature_pos | sigfeature_neg)


si = [ccpairs.strfsimilarity];

sum(isnan(si))

cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
markersize = 4;
linewidth = 2;
ymax = 0.15;

figure;
subplot(2,1,1);
hold on;

plot(si(sigfeature_pos), ccc_pos(sigfeature_pos), 'o', 'color', cmap(1,:), ...
   'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
   'markeredgecolor', cmap(1,:));

plot(si(sigfeature_neg), ccc_neg(sigfeature_neg), 'o', 'color', cmap(3,:), ...
   'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
   'markeredgecolor', cmap(3,:));


ytick = [0.0001 0.001 0.01 0.1 1];
xtick = [0:0.2:1];
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');
xlabel('STRF Similarity');
ylabel('Correlation');
legend('Facilitative', 'Suppressive');
tickpref;
axis([0 1 0.001 1]);



subplot(2,1,2);
data = si;
edges = 0:0.05:1;
n = histc(data, edges);
n = n ./ sum(n);
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
% center = edge2center(edges);
% n = n(1:end-1);
% hold on;
% plot(center, n, 'o-', 'color', cmap(j,:), ...
% 'markersize', markersize, 'markerfacecolor', cmap(j,:), ...
% 'markeredgecolor', cmap(j,:), 'linewidth', linewidth);
box off;
tickpref;
xtick = edges(1:2:end);
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = linspace(0,ymax,4);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([min(edges)-0.025 max(edges)+0.025]);
ylim([0 ymax]);
ylabel('Proportion');
xlabel('STRF Similarity');


fprintf('\n');
fprintf('STRF Similarity population\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', ...
   length(si), mean(si(~isnan(si))), ...
   std(si(~isnan(si) ) ) );
fprintf('MD = %.4f, MAD=%.4f, MX = %.4f\n', ...
   median(si(~isnan(si))), mad(si(~isnan(si))), ...
   max(si(~isnan(si))) );
fprintf('\n\n');


ccpairs_mn = mean(si);
ccpairs_sd = std(si);
ccpairs_se = ccpairs_sd / sqrt(length(si));




[r,p] = corrcoef(si(sigfeature_pos), log10(ccc_pos(sigfeature_pos)) );
fprintf('\nCCPairs: log10 CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigfeature_pos) )), log10(abs(ccc_pos(sigfeature_pos) )));
fprintf('CCPairs: log10 CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigfeature_pos), log10(ccc_pos(sigfeature_pos) ));
fprintf('CCPairs: log10 CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

% xfit = linspace(0.1, 1, 1000);
% [beta,s] = polyfit(ccpairs_si_sig, log10(ccpairs_rho_sig), 1); % linear regression
% yfit = polyval(beta,xfit);
% hp = plot(xfit, 10.^yfit, '-', 'color', cmap(1,:), 'linewidth', linewidth);

box off;
tickpref;
set(gcf,'position',[208   444   321   420]);
print_mfilename(mfilename);

   
return;







function plot_allpairs_spktype_strf_crosscorr_similarity(ccpairs)

ccpairs
pause


pd_pos = [ccpairs.pd_pos];
hw_pos = [ccpairs.hw_pos];
sigfeature_pos = [ccpairs.sigfeature_pos];
ccc_pos = [ccpairs.ccc_pos];


pd_neg = [ccpairs.pd_neg];
hw_neg = [ccpairs.hw_neg];
sigfeature_neg = [ccpairs.sigfeature_neg];
ccc_neg = [ccpairs.ccc_neg];





% Use new calculation of correlation coefficient if available
if ( isfield(ccpairs, 'ccc') )
   ccpairs_rho = [ccpairs.ccc];
else
   ccpairs_rho = [ccpairs.rho];
end
ccpairs_rho(ccpairs_rho<=0) = 0.0001;

ccpairs_si = [ccpairs.strfsimilarity];
ccpairs_sig = logical([ccpairs.significant]);

ccpairs_rho_sig = ccpairs_rho(ccpairs_sig);
ccpairs_si_sig = ccpairs_si(ccpairs_sig);

% Now plot the data according to cell type

spktype_data{1} = ccpairs_si;

spktype_label{1} = 'CCPairs';

cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
markersize = 4;
linewidth = 2;
ymax = 0.15;

figure;
subplot(2,1,1);
hold on;
index = find(ccpairs_rho > 0.0001);
ccpairs_rho = ccpairs_rho(index);
ccpairs_si = ccpairs_si(index);
plot(ccpairs_si, ccpairs_rho, 'ko', 'color', 'k', ...
   'markersize', markersize, 'markerfacecolor', 'k', ...
   'markeredgecolor', 'k');
% plot(ccpairs_si_sig, ccpairs_rho_sig, 'o', 'color', cmap(1,:), ...
%    'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
%    'markeredgecolor', cmap(1,:));
% 
ytick = [0.0001 0.001 0.01 0.1 1];
xtick = [0:0.2:1];
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');
xlabel('STRF Similarity');
ylabel('Correlation');
% set(gca,'xscale', 'log');
tickpref;
axis([0 1 0.0001 1]);



subplot(2,1,2);
data = spktype_data{1};
edges = 0:0.05:1;
n = histc(data, edges);
n = n ./ sum(n);
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
% center = edge2center(edges);
% n = n(1:end-1);
% hold on;
% plot(center, n, 'o-', 'color', cmap(j,:), ...
% 'markersize', markersize, 'markerfacecolor', cmap(j,:), ...
% 'markeredgecolor', cmap(j,:), 'linewidth', linewidth);
box off;
tickpref;
xtick = edges(1:2:end);
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = linspace(0,ymax,4);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([min(edges)-0.025 max(edges)+0.025]);
ylim([0 ymax]);
ylabel('Proportion');
xlabel('STRF Similarity');


fprintf('\n');
fprintf('STRF Similarity population\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', ...
   length(ccpairs_si), mean(ccpairs_si(~isnan(ccpairs_si))), ...
   std(ccpairs_si(~isnan(ccpairs_si) ) ) );
fprintf('MD = %.4f, MAD=%.4f, MX = %.4f\n', ...
   median(ccpairs_si(~isnan(ccpairs_si))), mad(ccpairs_si(~isnan(ccpairs_si))), ...
   max(ccpairs_si(~isnan(ccpairs_si))) );
fprintf('\n\n');


ccpairs_mn = mean(ccpairs_si);
ccpairs_sd = std(ccpairs_si);
ccpairs_se = ccpairs_sd / sqrt(length(ccpairs_si));




[r,p] = corrcoef(ccpairs_si_sig, log10(ccpairs_rho_sig));
fprintf('\nCCPairs: log10 CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs(ccpairs_si_sig)), log10(abs(ccpairs_rho_sig)));
fprintf('CCPairs: log10 CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(ccpairs_si, log10(ccpairs_rho));
fprintf('CCPairs: log10 CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

% xfit = linspace(0.1, 1, 1000);
% [beta,s] = polyfit(ccpairs_si_sig, log10(ccpairs_rho_sig), 1); % linear regression
% yfit = polyval(beta,xfit);
% hp = plot(xfit, 10.^yfit, '-', 'color', cmap(1,:), 'linewidth', linewidth);

box off;
tickpref;
set(gcf,'position',[208   444   321   420]);
print_mfilename(mfilename);

   
return;



function plot_allpairs_spktype_strf_crosscorr_cf_q_latency(ccpairs) 
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% -----------------------------------------------------------------------
%
%
%
% caa 4/21/11

cf1 = [ccpairs.cf1];
cf2 = [ccpairs.cf2];

q1 = [ccpairs.q1];
q2 = [ccpairs.q2];

latency1 = [ccpairs.latency1];
latency2 = [ccpairs.latency2];

% dur1 = [ccpairs.dur1];
% dur2 = [ccpairs.dur2];

ccpairs_rho = [ccpairs.ccc];
ccpairs_rho(ccpairs_rho<=0) = 0.0001;
ccpairs_sig = logical([ccpairs.significant]);

cf = [cf1(:) cf2(:)];
q = [q1(:) q2(:)];
latency = [latency1(:) latency2(:)];

% [ccpairs_cf, ccpairs_q, ccpairs_latency, ccpairs_dur1, ccpairs_dur2] = ...
%    sort_dur_ptparams(cf, q, latency, dur1, dur2);

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
edges = {linspace(0, 0.5, 11), linspace(0, 5, 11), linspace(0, 20, 11)};
tick = {[0.001 0.01 0.1 1 3 5 10 20], [0.5 1 2 4 8], [0:4:20]};
ticklabel = {[0.001 0.01 0.1 1 3 5 10 20], [0.5 1 2 4 8], [0:4:20]};
xlabel1 = {'BF (kHz)', 'Q', 'Latency (ms)'};
xlabel2 = {'BF Difference (oct)', 'Q Difference (oct)', 'Latency Difference (ms)'};

spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'cf', 'q', 'latency'};
celltype = {'ccpairs'};

ymaxdiff = [1 1 0.65];
ytick = {linspace(0,1,5), linspace(0,1,5), linspace(0,0.6,4)}; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 4;

figure;

for i = 1:length(datatype)

      dataspktype = eval(['ccpairs_' datatype{i}]);
 
      subplot(3,1,i);
      hold on;

      if ( i > 2 )
         dataspktype(:,1) = jitter( dataspktype(:,1) );
         dataspktype(:,2) = jitter( dataspktype(:,2) );
      end 

      [larger,smaller]= vs_largersmaller(dataspktype(:,1),dataspktype(:,2));

      plot(larger, smaller, ...
         'ko', 'markerfacecolor', 'k', 'markersize', 2);
      
      fprintf('\nDatatype: %s\n', datatype{i});
      plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');
      set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
      tickpref;
      set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
      set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
      xlim([xmin(i) xmax(i)]);
      ylim([ymin(i) ymax(i)]);
      xlabel(sprintf('%s %s', spktype_xlabel{1}, xlabel1{i}));
      ylabel(sprintf('%s %s', spktype_ylabel{1}, xlabel1{i}));
      

%       Didn't need absolute value
      if ( strcmp(datatype{i}, 'cf') )
         datadiff = abs( log2( dataspktype(:,1) ./ dataspktype(:,2) ) );
         [r, p] = corrcoef(log10(dataspktype(:,1)), log10(dataspktype(:,2)) );
         fprintf('%s  CF r = %.3f, p = %.4f\n', spktype_label{1}, r(1,2), p(1,2) );
      elseif ( strcmp(datatype{i}, 'latency') )
         datadiff = abs( dataspktype(:,1) - dataspktype(:,2) );
         [r, p] = corrcoef((dataspktype(:,1)), (dataspktype(:,2)) );
         fprintf('%s  Latency r = %.3f, p = %.4f\n', spktype_label{1}, r(1,2), p(1,2) );
      else % it must be 'q'
         datadiff = abs( log2( dataspktype(:,1) ./ dataspktype(:,2) ) );
         [r, p] = corrcoef(log10(dataspktype(:,1)), log10(dataspktype(:,2)) );
         fprintf('%s  Q r = %.3f, p = %.4f\n', spktype_label{1}, r(1,2), p(1,2) );
      end
   
      title(sprintf('N = %.0f, r = %.3f', size(dataspktype,1), r(2)));
%       subplot(2,2,4);
%       n = histc(datadiff, edges{i});
%       n = n ./ sum(n);
%       hold on;
%       plot(edges{i}, n, [marker{j} '-'], 'color', cmap(j,:), ...
%          'markersize', 4, 'markerfacecolor', cmap(j,:), ...
%          'markeredgecolor', cmap(j,:), 'linewidth', linewidth);
%       box off;
%       set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%       xtick = edges{i}(1:2:end);
%       set(gca,'xtick', xtick, 'xticklabel', xtick);
%       set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
%       range = max(edges{i}) - min(edges{i});
%       xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
%       ylim([0 ymaxdiff(i)]);
%       ylabel('Proportion');
%       xlabel(xlabel2{i});
% 
%       fprintf('\n%s pairs - %s:\n', spktype_label{j}, xlabel1{i});
%       fprintf('%s median = %.3f\n', spktype_xlabel{j}, median(dataspktype(:,1)));
%       fprintf('%s median = %.3f\n', spktype_ylabel{j}, median(dataspktype(:,2)));
%       fprintf('Signed Rank p = %.4f\n', signrank(dataspktype(:,1),dataspktype(:,2)));
%       fprintf('\n');

   set(gcf,'position',[20 60 320 750]);
   print_mfilename(mfilename);

end % for i

return;



function plot_allpairs_spktype_strf_crosscorr_strf_fr_pli(ccpairs)
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% -----------------------------------------------------------------------
%
%
%
% caa 4/21/11



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

sum(sigfeature_pos & sigfeature_neg)
sum(sigfeature_pos & ~sigfeature_neg)
sum(~sigfeature_pos & sigfeature_neg)
sum(sigfeature_pos | sigfeature_neg)



% Data
fr1 = [ccpairs.fr1];
fr2 = [ccpairs.fr2];

pli1 = [ccpairs.pli1];
pli2 = [ccpairs.pli2];

si1 = [ccpairs.sepindex1];
si2 = [ccpairs.sepindex2];

% ccpairs_rho = [ccpairs.rho];
ccpairs_rho = [ccpairs.ccc];
ccpairs_rho(ccpairs_rho<=0) = 0.0001;
ccpairs_sig = logical([ccpairs.significant]);

fr = [fr1(:) fr2(:)];
pli = [pli1(:) pli2(:)];
si = [si1(:) si2(:)];

ccpairs_fr = fr;
ccpairs_pli = pli;
ccpairs_si = si;

% [ccpairs_fr, ccpairs_pli, ccpairs_si, ccpairs_dur1, ccpairs_dur2] = ...
%    sort_dur_strfparams(fr, pli, si, dur1, dur2);

xmin = [0.1 0.01 0.2];
xmax = [30 1 1];
ymin = xmin;
ymax = xmax;
xscale = {'log', 'log', 'log'};
yscale = {'log', 'log', 'log'};
edges = {0:2.5:20, 0:0.05:0.4, 0:0.075:0.6};
edges = {linspace(0,100,21), linspace(0,0.4,17), 0:0.075:0.6};
tick = {[0.1 1 10 30], [0.01 0.1 1], [0:0.2:1]};
ticklabel = {[5 10 20 30], [0.5 1 2 4 8], [0:10:50]};
xlabel1 = {'Firing Rate', 'RPI', 'Separability Index'};
xlabel2 = {'FR Difference (Hz)', 'RPI Difference', 'Separability Index Diff'};


spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'fr', 'pli', 'si'};
celltype = {'ccpairs'};

dataxlabel = {'Firing Rate (Hz) Neuron 1', 'RPI Neuron 1' };
dataylabel = {'Firing Rate (Hz) Neuron 2', 'RPI Neuron 2' };

datadiffxlabel = {'FR Difference (Hz)','RPI Difference (Hz)'};
datadiffylabel = {'Cross Correlation Coefficient','Cross Correlation Coefficient' };

ymaxdiff = [0.7 0.41 0.65];
ytick = {0:0.1:0.7, linspace(0,0.4,5), linspace(0,0.6,4)}; 
xtick2 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 1;
yax = {[0.001 1], [0.001 1], [0.0001 2]};
xax = {[0.01 100],[0.0001 1] [0.0001 2]};

% x = {log10(logspace(log10(0.0001), log10(100), 1000)) };

x = {log10(logspace(log10(0.01), log10(100), 1000)), ...
log10(logspace(log10(0.0001), log10(1), 1000))};

subplotnum = [1 3 5];

figure;

for i = 1:2 %length(datatype)

   fignum = (i-1) + subplotnum;

   dataspktype = eval(['ccpairs_' datatype{i}]);
   rho = eval(['ccpairs_rho']);
   sig = eval(['ccpairs_sig']);
   [larger,smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
   datadiff = abs( larger - smaller );

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

     
   subplot(3,2,fignum(2));
   n = histc(datadiff, edges{i});
   n = n ./ sum(n);
   hb = bar(edges{i},n,'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
   hold on;
   box off;
   tickpref;
   xtick = edges{i}(1:2:end);
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   range = max(edges{i}) - min(edges{i});
   xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
   ylim([0 ymaxdiff(i)]);
   ylabel('Proportion');
   xlabel(xlabel2{i});

     
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

set(gcf,'position',[247   101   322   360]);
print_mfilename(mfilename);

% se = sd ./ sqrt(length(datadiff));

return;



function plot_allpairs_spktype_strf_crosscorr_fio_similarity(ccpairs)

% Use new calculation of correlation coefficient if available
if ( isfield(ccpairs, 'ccc') )
   ccpairs_rho = [ccpairs.ccc];
else
   ccpairs_rho = [ccpairs.rho];
end
ccpairs_rho(ccpairs_rho<=0) = 0.0001;

ccpairs_si = [ccpairs.fiosimilarity];
ccpairs_sig = logical([ccpairs.significant]);

ccpairs_rho_sig = ccpairs_rho(ccpairs_sig);
ccpairs_si_sig = ccpairs_si(ccpairs_sig);


% Now plot the data according to cell type

spktype_data{1} = ccpairs_si;

spktype_label{1} = 'CCPairs';

cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
markersize = 4;
linewidth = 2;

figure;
for j = 1:length(spktype_label)
   data = spktype_data{j};
   edges = 0.4:0.025:1;
   n = histc(data, edges);
   n = n ./ sum(n);
   hold on;
   hb = bar(edges, n, 'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
end % (for j)

box off;
tickpref;
xtick = edges(1:2:end);
set(gca,'xtick', xtick, 'xticklabel', xtick);
ytick = 0:0.1:0.5;
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([min(edges) max(edges)]);
ylim([0 0.5]);
ylabel('Proportion');

ccpairs_mn = mean(ccpairs_si);
ccpairs_md = median(ccpairs_si);
ccpairs_sd = std(ccpairs_si);
ccpairs_se = ccpairs_sd / sqrt(length(ccpairs_si));

title(sprintf('Nonlinearity Similarity Indices\nMn=%.2f,Md=%.2f,SD=%.2f,SE=%.2f', ...
ccpairs_mn, ccpairs_md, ccpairs_sd, ccpairs_se))

fprintf('Nonlinearity Similarity Indices\nMn=%.2f,Md=%.2f,SD=%.2f,SE=%.2f', ...
ccpairs_mn, ccpairs_md, ccpairs_sd, ccpairs_se);

fprintf('\n');

box off;
tickpref;

set(gcf,'position', [1359 745 321 233]);
print_mfilename(mfilename);

return;



function plot_allpairs_spktype_strf_crosscorr_fio_asi_similarity(ccpairs)

asi1 = [ccpairs.fioasi1];
asi2 = [ccpairs.fioasi2];
asi = [asi1(:) asi2(:)];
ccpairs_asi = asi;
ccpairs_si = [ccpairs.fiosimilarity];

% Use new calculation of correlation coefficient if available
if ( isfield(ccpairs, 'ccc') )
   ccpairs_rho = [ccpairs.ccc];
else
   ccpairs_rho = [ccpairs.rho];
end
ccpairs_rho(ccpairs_rho<=0) = 0.0001;
ccpairs_sig = logical([ccpairs.significant]);




xmin = [0.5];
xmax = [1];
ymin = [0.5];
ymax = [1];
xscale = 'linear';
yscale = 'linear';
edges = 0:0.025:0.5;
tick = [0:0.1:1];
ticklabel = [0:0.1:1];
xlabel1 = {'Asymmetry Index', 'Asymmetry Index', 'Asymmetry Index'};
xlabel2 = {'ASI Difference', 'ASI Difference', 'ASI Difference'};
spktype_label = 'CCPairs';
spktype_xlabel = 'Neuron 1';
spktype_ylabel = 'Neuron 2';


ymaxdiff = [0.33];
ytick = linspace(0,1,11); 
ytick2 = [0.01 0.1 1]; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 1;
yax = [0.001 1];
xax = [0.001 1];
x = {linspace(0,1,1000)};
x = log10(logspace(log10(0.0001), log10(1), 1000));

figure;

dataspktype = ccpairs_asi; %eval(['ccpairs_' datatype{i}]);
[larger,smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
datadiff = abs( larger - smaller );


subplot(1,3,1);
hold on;
plot(larger, smaller, 'ko', 'markerfacecolor', 'k', 'markersize', markersize);
plot([xmin xmax], [ymin ymax], 'k-');
set(gca,'xscale', xscale, 'yscale', yscale);
tickpref;
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('ASI Neuron 1');
ylabel('ASI Neuron 2');


fprintf('\n');
fprintf('\nAsymmetry Index\n');
[r,p] = corrcoef( (larger), (smaller) );
fprintf('r=%.3f, p = %.3f\n', r(2),p(2));




subplot(1,3,2);
n = histc(datadiff, edges);
n = n ./ sum(n);
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', 0:0.1:0.5, 'xticklabel', 0:0.1:0.5);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([0 0.5]);
ylim([0 ymaxdiff]);
ylabel('Proportion');
xlabel('ASI Difference');


[pval, data_median] = pairs_permutation_test(larger, smaller);
fprintf('\nAsymmetry Index Difference\n');
fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
   length(datadiff),mean(datadiff), std(datadiff), std(datadiff)./sqrt(length(datadiff)));
fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
   median(datadiff), mad(datadiff), max(datadiff));
fprintf('ASI Permutation pval: %.4f\n', pval);
fprintf('\n\n');




subplot(1,3,3);
edges = 0.4:0.025:1;
n = histc(ccpairs_si, edges);
n = n ./ sum(n);
hold on;
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
box off;
tickpref;
xtick = edges(1:2:end);
set(gca,'xtick', 0.4:0.1:1, 'xticklabel', 0.4:0.1:1);
ytick = 0:0.1:0.5;
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([min(edges) max(edges)]);
ylim([0 0.5]);
ylabel('Proportion');
xlabel('Nonlinearity Similarity');
box off;
tickpref;

fprintf('\nFIO Similarity\n');
fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
   length(ccpairs_si),mean(ccpairs_si), std(ccpairs_si), std(ccpairs_si)./sqrt(length(ccpairs_si)));
fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
   median(ccpairs_si), mad(ccpairs_si), max(ccpairs_si));


set(gcf,'position', [400 800 438 118]);
print_mfilename(mfilename);



return;



function plot_allpairs_spktype_strf_crosscorr_fio_asi(ccpairs)


asi1 = [ccpairs.fioasi1];
asi2 = [ccpairs.fioasi2];

% dur1 = [ccpairs.dur1];
% dur2 = [ccpairs.dur2];

% Use new calculation of correlation coefficient if available
if ( isfield(ccpairs, 'ccc') )
   ccpairs_rho = [ccpairs.ccc];
else
   ccpairs_rho = [ccpairs.rho];
end
ccpairs_rho(ccpairs_rho<=0) = 0.0001;
ccpairs_sig = logical([ccpairs.significant]);

asi = [asi1(:) asi2(:)];
ccpairs_asi = asi;
% [ccpairs_asi, ccpairs_dur1, ccpairs_dur2] = sort_dur_asi(asi, dur1, dur2);

xmin = [0.5 0.4 0.4];
xmax = [1 1 1];
ymin = [0.5 0.4 0.4];
ymax = [1 1 1];
xscale = {'linear', 'linear', 'linear'};
yscale = {'linear', 'linear', 'linear'};
edges = {0:0.025:0.5};
tick = {[0:0.1:1], [0:0.1:1], [0:0.1:1]};
ticklabel = {[0:0.1:1], [0:0.1:1], [0:0.1:1]};
xlabel1 = {'Asymmetry Index', 'Asymmetry Index', 'Asymmetry Index'};
xlabel2 = {'ASI Difference', 'ASI Difference', 'ASI Difference'};


spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'asi', 'asi', 'asi'};
celltype = {'ccpairs'};

ymaxdiff = [0.33 0.8 0.65];
ytick = {linspace(0,1,11), linspace(0,0.8,5), linspace(0,0.6,4)}; 
ytick2 = {[0.01 0.1 1]}; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 2;
yax = {[0.001 1], [0.001 1]};
xax = {[0.001 1], [0.001 1]};

% x = {log10(logspace(log10(0.0001), log10(100), 1000)) };

x = {linspace(0,1,1000)};

x = {log10(logspace(log10(0.0001), log10(1), 1000))};


for i = 1:1 %2%length(datatype)

   figure;

   for j = 1:length(celltype)

      dataspktype = eval([celltype{j} '_' datatype{i}]);
      datadiff = abs( dataspktype(:,1) - dataspktype(:,2) );

      rho = eval([celltype{j} '_rho']);
      sig = eval([celltype{j} '_sig']);


      dataspktype = eval(['ccpairs_' datatype{i}]);
      rho = eval(['ccpairs_rho']);
      sig = eval(['ccpairs_sig']);

      [larger,smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));

      subplot(2,1,1);
      hold on;
      plot(larger, smaller, ...
         'ko', 'markerfacecolor', 'k', 'markersize', 1.5);

      plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');
      set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
      tickpref;
      set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
      set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
      xlim([xmin(i) xmax(i)]);
      ylim([ymin(i) ymax(i)]);
      xlabel(sprintf('%s %s', spktype_xlabel{j}, xlabel1{i}));
      ylabel(sprintf('%s %s', spktype_ylabel{j}, xlabel1{i}));
      title(sprintf('N = %.0f', size(dataspktype,1)));


      subplot(2,1,2);
      n = histc(datadiff, edges{i});
      n = n ./ sum(n);
      hb = bar(edges{i}, n, 'histc');
      set(hb, 'facecolor', [0.7 0.7 0.7]);
      asi_mean_diff = mean(datadiff);
      asi_median_diff = median(datadiff);
      title(sprintf('ASI Diff, Mean = %.2f, Median = %.2f',asi_mean_diff,asi_median_diff))

      box off;
      set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
      xtick = edges{i}(1:2:end);
      set(gca,'xtick', xtick, 'xticklabel', xtick);
      set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
      range = max(edges{i}) - min(edges{i});
      xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
      ylim([0 ymaxdiff(i)]);
      ylabel('Proportion');
      xlabel(xlabel2{i});



%       subplot(2,4,j+4);
%       hold on;
%       plot(datadiff, rho, 'o', 'color', cmap(j,:), ...
%          'markersize', markersize, 'markerfacecolor', cmap(j,:), ...
%          'markeredgecolor', cmap(j,:));
%       set(gca,'xtick', ytick2{i}, 'xticklabel', ytick2{i});
%       set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
% 
%       tickpref;
%       box off;
%       xlim(xax{i});
%       ylim(yax{i});
% 
%       xlabel(xlabel2{i});
%       ylabel('Correlation Coefficient');
%       [r,p] = corrcoef(log10(datadiff(:)), log10(rho(:)));
%       title(sprintf('%s, r = %.3f, p = %.3f', ...
%          spktype_label{j}, r(1,2), p(1,2)));
% 
% 
%       [beta,s] = polyfit(log10(datadiff(:)), log10(rho(:)), 1); % linear regression
%       xfit = x{i};
%       yfit = polyval(beta,xfit);
%       hp = plot(10.^xfit, 10.^yfit, '-', 'color', cmap(j,:), 'linewidth', linewidth);
% 
%       set(gca,'yscale', 'log');
%       set(gca,'xscale', 'log');
% 
%       fprintf('\n%s pairs - %s:\n', spktype_label{j}, xlabel1{i});
%       fprintf('%s median = %.3f, mad = %.3f\n', ...
%          spktype_xlabel{j}, median(dataspktype(:,1)), mad(dataspktype(:,1)) );
%       fprintf('%s mean = %.3f, SD = %.3f\n', ...
%          spktype_xlabel{j}, mean(dataspktype(:,1)), std(dataspktype(:,1)) );
% 
%       fprintf('%s median = %.3f, mad = %.3f\n', ...
%          spktype_ylabel{j}, median(dataspktype(:,2)), mad(dataspktype(:,2)));
%       fprintf('%s mean = %.3f, SD = %.3f\n', ...
%          spktype_ylabel{j}, mean(dataspktype(:,2)), std(dataspktype(:,2)) );
% 
%       fprintf('Signed Rank p = %.4f\n', ranksum(dataspktype(:,1),dataspktype(:,2)));
%       fprintf('\n');

   end % (for j)

%    legend(spktype_label{1:3});
   set(gcf,'position', [150   250   321 470]);
   print_mfilename(mfilename);

end % for i

return;










function plot_allpairs_spktype_strf_crosscorr_mtf_btmf(ccpairs)

% Pairs Data
btmf1 = [ccpairs.btmf1];
btmf2 = [ccpairs.btmf2];

bsmf1 = [ccpairs.bsmf1];
bsmf2 = [ccpairs.bsmf2];

% dur1 = [ccpairs.dur1];
% dur2 = [ccpairs.dur2];


% Use new calculation of correlation coefficient if available
if ( isfield(ccpairs, 'ccc') )
   ccpairs_rho = [ccpairs.ccc];
else
   ccpairs_rho = [ccpairs.rho];
end
ccpairs_rho(ccpairs_rho<=0) = 0.0001;
ccpairs_sig = logical([ccpairs.significant]);

btmf = [btmf1(:) btmf2(:)];
bsmf = [bsmf1(:) bsmf2(:)];

ccpairs_btmf = btmf;
ccpairs_bsmf = bsmf;

% [ccpairs_btmf, ccpairs_dur1, ccpairs_dur2] = sort_dur_asi(btmf, dur1, dur2);
% [ccpairs_bsmf, ccpairs_dur1, ccpairs_dur2] = sort_dur_asi(bsmf, dur1, dur2);

spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';
datatype = {'btmf', 'bsmf'};
celltype = {'ccpairs'};

xmin = [12.5 0.05];
xmax = [400 2.1];
ymin = xmin;
ymax = xmax;
xscale = {'log', 'log', 'log'};
yscale = {'log', 'log', 'log'};
edges = {linspace(0, 300, 21), linspace(0, 2, 21)};
tick = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
ticklabel = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
xlabel1 = {'bTMF (Hz)', 'bSMF (cyc/oct)'};
xlabel2 = {'bTMF Difference (Hz)', 'bSMF Difference (cyc/oct)'};

ymaxdiff = [.27 0.5];
ytick = {0:0.05:0.25, 0:0.1:0.5, linspace(0,0.6,4)}; 
xtick2 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1 2]}; 
ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1]}; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o', 'o', 'o'};
linewidth = 2;
markersize = 2;
yax = {[0.001 1], [0.001 1]};
xax = {[0.01 400],[0.001 2]};

% x = {log10(logspace(log10(0.0001), log10(100), 1000)) };

x = {log10(logspace(log10(0.01), log10(100), 1000)), ...
log10(logspace(log10(0.001), log10(10), 1000))};

for i = 1:length(datatype)

   figure;

   dataspktype = eval(['ccpairs_' datatype{i}]);
   rho = eval(['ccpairs_rho']);
   sig = eval(['ccpairs_sig']);
   [larger, smaller] = vs_largersmaller(dataspktype(:,1),dataspktype(:,2));
   datadiff = abs( dataspktype(:,1) - dataspktype(:,2) );

   subplot(3,1,1);
   hold on;
   plot(larger, smaller, 'ko', 'markerfacecolor', 'k', 'markersize', 2);
   plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');
   [r,p] = corrcoef(log10(dataspktype(:,1)), log10(dataspktype(:,2)));
   set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
   tickpref;
   set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
   set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
   xlim([xmin(i) xmax(i)]);
   ylim([ymin(i) ymax(i)]);
   xlabel(sprintf('%s %s', spktype_xlabel{1}, xlabel1{i}));
   ylabel(sprintf('%s %s', spktype_ylabel{1}, xlabel1{i}));
   title(sprintf('N = %.0f, r = %.3f, p = %.3f', size(dataspktype,1),...
      r(2), p(2)));


   subplot(3,1,2);
   n = histc(datadiff, edges{i});
   n = n ./ sum(n);
   hb = bar(edges{i}, n, 'histc');
   set(hb, 'facecolor', [0.7 0.7 0.7]);
   box off;
   tickpref; %set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   xtick = edges{i}(1:2:end);
   %       set(gca,'xtick', xtick, 'xticklabel', xtick);
   set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
   range = max(edges{i}) - min(edges{i});
   %       xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
   xlim([min(edges{i}) max(edges{i})]);
   ylim([0 ymaxdiff(i)]);
   ylabel('Proportion');
   xlabel(xlabel2{i});

   subplot(3,1,3);
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
   [r,p] = corrcoef((datadiff(:)), log10(rho(:)));
   title(sprintf('%s, r = %.3f, p = %.3f', spktype_label{1}, r(1,2), p(1,2)));
   set(gca,'yscale', 'log');
   set(gca,'xscale', 'log');

   set(gcf,'position',[50    60   321   750]);
   print_mfilename(mfilename);

end % for i



return




function [asi, dur1, dur2] = sort_dur_asi(asi, dur1, dur2)
%sort_dur_strfparams
%

for i = 1:length(dur1)

   if ( dur1(i) > dur2(i) )

      a = dur1(i);
      b = dur2(i);
      dur1(i) = b;
      dur2(i) = a;

      a = asi(i,1);
      b = asi(i,2);
      asi(i,1) = b;
      asi(i,2) = a;

   end

end % (for i)

return;




function [cf, q, latency, dur1, dur2] = sort_dur_ptparams(cf, q, latency, dur1, dur2)
% Now sort the data so that dur1 has the lowest durations, dur2 the
% longest. Swap the data values accordingly

for i = 1:length(dur1)

   if ( dur1(i) > dur2(i) )

      a = dur1(i);
      b = dur2(i);
      dur1(i) = b;
      dur2(i) = a;

      a = cf(i,1);
      b = cf(i,2);
      cf(i,1) = b;
      cf(i,2) = a;
   
      if ( max(cf(:)) > 1000 )
         cf = cf ./ 1000;
      end

      a = q(i,1);
      b = q(i,2);
      q(i,1) = b;
      q(i,2) = a;

      a = latency(i,1);
      b = latency(i,2);
      latency(i,1) = b;
      latency(i,2) = a;

   end

end % (for i)

return;




function [fr, pli, si, dur1, dur2] = sort_dur_strfparams(fr, pli, si, dur1, dur2)
%sort_dur_strfparams
%

for i = 1:length(dur1)

   if ( dur1(i) > dur2(i) )

      a = dur1(i);
      b = dur2(i);
      dur1(i) = b;
      dur2(i) = a;

      a = fr(i,1);
      b = fr(i,2);
      fr(i,1) = b;
      fr(i,2) = a;

      a = pli(i,1);
      b = pli(i,2);
      pli(i,1) = b;
      pli(i,2) = a;

      a = si(i,1);
      b = si(i,2);
      si(i,1) = b;
      si(i,2) = a;

   end

end % (for i)

return;






