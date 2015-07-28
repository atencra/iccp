function [rbarpop, dtpop, infopop, rbarpairs, dtpairs, infopairs] = ...
   vs_iccpairs_plot_info_calc_pop(rbarpop, dtpop, infopop, rbarpairs, dtpairs, infopairs)
% NEEDS TO BE MODULARIZED
%
% info_calculation: Calculates information values and nmse values from an
% STRF response file. (saved as resp_final)
%
% [info_nmse, sigma_nmse, nmse] = info_calculation(resp, graphs, popgraphs)
%
% Inputs:
% -----------------------------
% 1) resp: STRF response file, saved as -*final-resp. 
% 2) graphs: boolean input. 1 if graphs are desired, 0 if not.
% 3) popgraphs: boolean input. 1 if population graphs are desired, 0 if
%       not.
%
% Outputs:
%   ----------------------------
%   IF GRAPHS is 1 (plots are for each individual neuron):
%
%   1) Generates a plot of Information vs. all the dt values.
%   2) Generates a PSTH, raster plot, and image of the stimulus 
%             spectrogram, using the value of dt where the I vs dt graph 
%             becomes nonlinear.
%   3) Plots all of its information values, divided into intervals of 
%             width [[0.1 0.5 1 2 3 4 6 12], against log(1/those T values).
%   4) Plots sigma against 1/sqrt(those T values).
%
%   IF POPGRAPHS is 1 (plots are for all the neurons)
%     1) Plots the population histogram of normalized mean square error, 
%             information, and picked dt values.
%
%  Finally, performs pairwise analysis on information, sigma values, and
%  generates scatterplots of those.

% Things to plot:
% 1. Raster
% 2. PSTH
% 3. Information vs. time bin size
% 4. Information vs. data set size for multiple estimates
% 5. SD vs. data set size


% If there's no inputs, then extract the data from files
if ( nargin ~= 6 )
   [rbarpop, dtpop, infopop, rbarpairs, dtpairs, infopairs] = ...
      vs_iccpairs_get_info_calc_pop_data;
end % if ( nargin )


close all;

% Population figures
figure;

% Population histograms
subplot(3,2,1);
vs_iccpairs_plot_rbar_pop(rbarpop);

subplot(3,2,3);
vs_iccpairs_plot_dt_pop(dtpop);

subplot(3,2,5);
vs_iccpairs_plot_info_pop(infopop);


% Population comparisons
subplot(3,2,2);
vs_iccpairs_plot_rbar_dt_pop(rbarpop, dtpop);

subplot(3,2,4);
vs_iccpairs_plot_rbar_info_pop(rbarpop, infopop);

subplot(3,2,6);
vs_iccpairs_plot_dt_info_pop(dtpop, infopop);

% set(gcf,'position', [1070 118 438 539]);
print_mfilename(mfilename);


close all;

% Pairwise figures
figure;

subplot(3,1,1);
vs_iccpairs_plot_rbar_pairs(rbarpairs);

subplot(3,1,2);
vs_iccpairs_plot_dt_pairs(dtpairs);

subplot(3,1,3);
vs_iccpairs_plot_info_pairs(infopairs);

set(gcf,'position', [577 171 321 759]);
print_mfilename(mfilename);

% figure;
% vs_iccpairs_plot_rbar_diff_pairs(rbarpairs);

% figure;
vs_iccpairs_plot_dt_diff_pairs(dtpairs);

figure;
vs_iccpairs_plot_info_diff_pairs(infopairs);



return;







function [rbarpop, dtpop, infopop, rbarpairs, dtpairs, infopairs] = ...
   vs_iccpairs_get_info_calc_pop_data

   d = dir('*-respcmbcmb-info.mat');

   infopairs = [];
   dtpairs = [];
   rbarpairs = [];

   infopop = [];
   dtpop = [];
   rbarpop = [];

   for n = 1:length(d)

      filename = d(n).name;
      s = load(filename, 'iresp');
      iresp = s.iresp;

      % Get the data individually
      for i = 1:length(iresp)
         if ( ~isempty(iresp(i).psth) )
            infopop = [infopop; iresp(i).ifracdurmn(end)];
            dtpop = [dtpop; iresp(i).dtoptim];
            rbarpop = [rbarpop; iresp(i).rbar];
         else
            infopop = [infopop; NaN];
            dtpop = [dtpop; NaN];
            rbarpop = [rbarpop; NaN];
         end
      end % (for i)


      % Get the data pair-wise
      chan = [iresp.chan];
      chan_unique = unique(chan);

      for i = 1:length(chan_unique)

         index = find( chan_unique(i) == chan );

         if ( length(index)>1 )

            cmb = nchoosek(index, 2); % determine all possible pairwise combinations

            [nr, nc] = size(cmb);

            for j = 1:nr

               index1 = cmb(j,1);
               index2 = cmb(j,2);

               iresp1 = iresp(index1);
               iresp2 = iresp(index2);

               % First save the statistics for later

               if ( ~isempty(iresp1.psth) && ~isempty(iresp2.psth)  )
                  infopairs = [infopairs; iresp1.ifracdurmn(end) iresp2.ifracdurmn(end)];
                  dtpairs = [dtpairs; iresp1.dtoptim iresp2.dtoptim];
                  rbarpairs = [rbarpairs; iresp1.rbar iresp2.rbar];

               elseif ( isempty(iresp1.psth) && ~isempty(iresp2.psth)  )
                  infopairs = [infopairs; NaN iresp2.ifracdurmn(end)];
                  dtpairs = [dtpairs; NaN iresp2.dtoptim];
                  rbarpairs = [rbarpairs; NaN iresp2.rbar];

               elseif ( ~isempty(iresp1.psth) && isempty(iresp2.psth)  )
                  infopairs = [infopairs; iresp1.ifracdurmn(end) NaN];
                  dtpairs = [dtpairs; iresp1.dtoptim NaN];
                  rbarpairs = [rbarpairs; iresp1.rbar NaN];

               elseif ( isempty(iresp1.psth) && isempty(iresp2.psth)  )
                  infopairs = [infopairs; NaN NaN];
                  dtpairs = [dtpairs; NaN NaN];
                  rbarpairs = [rbarpairs; NaN NaN];
               end

            end % (for j)

         end % (if/else)

       end % (for i)

   end % ( for n = 1:length(d) )

return;





function vs_iccpairs_plot_rbar_pop(rbarpop)

hold on;
edges = linspace(0,60,16);
edges = 0:2.5:60;
n = histc(rbarpop, edges);
n = n ./ sum(n);
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
xlim([0 max(edges)])
ylim([0 1.05*max(n)]);
box off;
tickpref;
ylabel('Proportion');
xlabel('Ave Firing Rate (sp/s)');
xtick = [0:10:60];
set(gca,'xtick', xtick, 'xticklabel', xtick);
text(min(xlim), 1.1*max(ylim), 'A');

fprintf('\n');
fprintf('Rbar population\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', length(rbarpop), mean(rbarpop(~isnan(rbarpop))), std(rbarpop(~isnan(rbarpop) ) ) );
fprintf('MD = %.4f, MAD=%.4f, MX = %.4f\n', median(rbarpop(~isnan(rbarpop))), mad(rbarpop(~isnan(rbarpop))), max(rbarpop(~isnan(rbarpop))) );
fprintf('\n\n');


return;



function vs_iccpairs_plot_dt_pop(dtpop)

hold on;
edges = linspace(0,0.030,16);
edges = 0:0.002:0.024;
edges = edges * 1000;
n = histc(1000*dtpop, edges);
n = n ./ sum(n);
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
xlim([0 max(edges)])
ylim([0 1.05*max(n)]);
box off;
tickpref;
ylabel('Proportion');
xlabel('Bin Size (ms)');
xtick = [0:5:20];
set(gca,'xtick', xtick, 'xticklabel', xtick);
text(min(xlim), 1.1*max(ylim), 'B');

dtpop = dtpop * 1000;
fprintf('\n');
fprintf('DT population\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', length(dtpop), mean(dtpop(~isnan(dtpop))), std(dtpop(~isnan(dtpop) ) ) );
fprintf('MD = %.4f, MAD=%.4f, MX = %.4f\n', median(dtpop(~isnan(dtpop))), mad(dtpop(~isnan(dtpop))), max(dtpop(~isnan(dtpop))) );
fprintf('\n\n');


return;



function vs_iccpairs_plot_info_pop(infopop)

hold on;
edges = linspace(0,10,21);
n = histc(infopop, edges);
n = n ./ sum(n);
hb = bar(edges, n, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
xlim([0 max(edges)])
ylim([0 1.05*max(n)]);
box off;
tickpref;
ylabel('Proportion');
xlabel('Information (bits)');
xtick = [0:2.5:10];
set(gca,'xtick', xtick, 'xticklabel', xtick);
text(min(xlim), 1.1*max(ylim), 'C');

fprintf('\n');
fprintf('Info population\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', length(infopop), mean(infopop(~isnan(infopop))), std(infopop(~isnan(infopop) ) ) );
fprintf('MD = %.4f, MAD=%.4f, MX = %.4f\n', median(infopop(~isnan(infopop))), mad(infopop(~isnan(infopop))), max(infopop(~isnan(infopop))) );
fprintf('\n\n');

return;






function vs_iccpairs_plot_rbar_dt_pop(rbarpop, dtpop)

hold on;
dtpop = 1000 * dtpop;
plot(rbarpop, dtpop, 'ko', 'markersize', 3);
ylim([0 max(dtpop)])
xlim([0 max(rbarpop)]);
box off;
tickpref;
xlabel('Ave Firing Rate (sp/s)');
ylabel('Bin Size (ms)');
set(gca,'xscale', 'log');

xtick = [0.10 1 10 50];
xlim([min(xtick) max(rbarpop)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);

ytick = 1000*[0 0.005 0.01 0.015 0.02 0.025];
ylim([min(ytick) max(ytick)]);
set(gca,'ytick', ytick, 'yticklabel', ytick);

text(min(xlim), 1.1*max(ylim), 'D');


index = ~isnan(rbarpop) & ~isnan(dtpop);
fprintf('\n');
fprintf('Rbar vs. DT population\n');
[r,p] = corrcoef(log10(rbarpop(index)), dtpop(index));
fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
fprintf('\n\n');


return;



function vs_iccpairs_plot_rbar_info_pop(rbarpop, infopop)

hold on;

plot(rbarpop, infopop, 'ko', 'markersize', 3);
ylim([0 max(infopop)])
box off;
tickpref;
xlabel('Ave Firing Rate (sp/s)');
ylabel('Information (bits)');

xtick = [0.10 1 10 50];
xlim([min(xtick) max(rbarpop)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'xscale', 'log');

ytick = [0.2 1 10];
ylim([min(ytick) max(ytick)]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');

text(min(xlim), 1.1*max(ylim), 'E');

index = ~isnan(rbarpop) & ~isnan(infopop);
fprintf('\n');
fprintf('Rbar vs. Info population\n');
[r,p] = corrcoef(log10(rbarpop(index)), log10(infopop(index)));
fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
fprintf('\n\n');


return;



function vs_iccpairs_plot_dt_info_pop(dtpop, infopop)

hold on;
dtpop = dtpop * 1000;
plot(dtpop, infopop, 'ko', 'markersize', 3);
xlim([0 max(dtpop)])
ylim([0 max(infopop)]);
box off;
tickpref;
ylabel('Information (bits/sp)');
xlabel('Bin Size (ms)');

xtick = 1000*[0 0.005 0.01 0.015 0.02 0.025];
xlim([min(xtick) max(xtick)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);

ytick = [0.2 1 10];
ylim([min(ytick) max(ytick)]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');

text(min(xlim), 1.1*max(ylim), 'F');


index = ~isnan(infopop) & ~isnan(dtpop);
fprintf('\n');
fprintf('Info vs. DT population\n');
[r,p] = corrcoef(log10(infopop(index)), dtpop(index));
fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
fprintf('\n\n');

return;








function vs_iccpairs_plot_rbar_diff_pairs(rbarpairs)


xmin = [0.1 0.01 0.2];
xmax = [30 1 1];
ymin = xmin;
ymax = xmax;
xscale = {'log'};
yscale = {'log'};
edges = {0:2.5:20};
edges = {linspace(0,80,21)};
tick = {[0.1 1 10 30]};
ticklabel = {[5 10 20 30]};
xlabel1 = {'Firing Rate'};
xlabel2 = {'FR Difference (Hz)'};


spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';

datatype = {'fr'};
celltype = {'ccpairs'};

dataxlabel = {'Firing Rate (Hz) Neuron 1'};
dataylabel = {'Firing Rate (Hz) Neuron 2'};

datadiffxlabel = {'FR Difference (Hz)'};
datadiffylabel = {'Cross Correlation Coefficient'};

ymaxdiff = [0.8 0.41 0.65];
ytick = {0:0.1:0.8}; 
xtick2 = {[0.01 0.1 1 10 100]}; 
ytick2 = {[0.001 0.01 0.1 1]}; 
cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
marker = {'o'};
linewidth = 2;
markersize = 1;
yax = {[0.001 1]};
xax = {[0.01 100]};

[larger,smaller] = vs_largersmaller(rbarpairs(:,1), rbarpairs(:,2));
datadiff = abs( larger - smaller );

data_median = nanmedian(datadiff);
data_mean = nanmean(datadiff);


subplot(1,2,1);
n = histc(datadiff, edges{1});
n = n ./ sum(n);
hb = bar(edges{1},n,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
hold on;
box off;
tickpref;
xtick = edges{1}(1:2:end);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick{1}, 'yticklabel', ytick{1});
range = max(edges{1}) - min(edges{1});
xlim([min(edges{1})-0.025*range max(edges{1})+0.025*range]);
ylim([0 ymaxdiff(1)]);
ylabel('Proportion');
xlabel(xlabel2{1});


% Now calculate what we would expect if the values were randomly assigned.
nreps = 1000;
rand_median = zeros(1,nreps);
rand_mean = zeros(1,nreps);

for i = 1:nreps
   index1 = randperm(length(larger));
   index2 = randperm(length(larger));
   temp1 = larger(index1);
   temp2 = smaller(index2);
   datadiff = abs(temp1 - temp2);
   rand_median(i) = nanmedian(datadiff);
   rand_mean(i) = nanmean(datadiff);
end % (for i)

index = find(rand_median < data_median);
pval = length(index) / nreps;


subplot(1,2,2);
hold on;
[n,x] = hist(rand_median,50);
edges = center2edge(x);
n = histc(rand_median, edges);
n = n ./ sum(n);
hb = bar(edges,n,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
set(hb, 'edgecolor', 'none');
tickpref;
minmin = min([data_median min(rand_median)]);
maxmax = max([data_median max(rand_median)]);
range = maxmax - minmin;
xlim([minmin-0.05*range maxmax+0.05*range]);
ylim([0 max(n)]);
plot([data_median data_median], ylim, 'k-', 'linewidth', 2);
xlabel('Median FR Difference');
title(sprintf('pval = %.4f, n = %.0f', pval, nreps));
set(gcf,'position', [577 171 438 169]);
print_mfilename(mfilename);



return;



function vs_iccpairs_plot_dt_diff_pairs(dtpairs)


dtpairs = 1000 * dtpairs; % convert from seconds to ms

edges = linspace(0,20,21);
xlabel1 = {'Bin Size'};
xlabel2 = {'Bin Size Difference (ms)'};

datadiffxlabel = {'Bin Size Difference (ms)'};
ymaxdiff = [0.41 0.65];
ytick1 = 0:0.1:1; 
xtick2 = {[0.1 1 10 100]}; 
ytick2 = {[0.01 0.1 1]}; 

[larger,smaller] = vs_largersmaller(dtpairs(:,1), dtpairs(:,2));
datadiff = abs( larger - smaller );
data_median = nanmedian(datadiff);
data_mean = nanmean(datadiff);


subplot(1,2,1);
n = histc(datadiff, edges);
n = n ./ sum(n);
hb = bar(edges,n,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
hold on;
box off;
tickpref;
xtick = edges(1:2:end);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick1, 'yticklabel', ytick1);
range = max(edges) - min(edges);
xlim([min(edges)-0.025*range max(edges)+0.025*range]);
ylim([0 1.05*max(n)]);
ylabel('Proportion');
xlabel('Bin Size Difference (ms)');


% Now calculate what we would expect if the values were randomly assigned.
nreps = 1000;
rand_median = zeros(1,nreps);
rand_mean = zeros(1,nreps);

for i = 1:nreps
   index1 = randperm(length(larger));
   index2 = randperm(length(larger));
   temp1 = larger(index1);
   temp2 = smaller(index2);
   datadiff = abs(temp1 - temp2);
   rand_median(i) = nanmedian(datadiff);
   rand_mean(i) = nanmean(datadiff);
end % (for i)
rand_median;

index = find(rand_median < data_median);
pval = length(index) / nreps;


% f1 = figure;
% hist(rand_median,25);
% min(rand_median)
% max(rand_median)
% pause
% close(f1);

subplot(1,2,2);
hold on;
[n,x] = hist(rand_median,50);
edges = center2edge(x);
edges;
n = histc(rand_median, edges);
n = n ./ sum(n);
hb = bar(edges,n,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
set(hb, 'edgecolor', 'none');
tickpref;
minmin = min([data_median min(rand_median)]);
maxmax = max([data_median max(rand_median)]);
range = maxmax - minmin;
xlim([minmin-0.05*range maxmax+0.05*range]);
ylim([0 max(n)]);
plot([data_median data_median], ylim, 'k-', 'linewidth', 2);
xlabel('Median Bin Size Difference');
title(sprintf('pval = %.4f, n = %.0f', pval, nreps));
set(gcf,'position', [577 171 438 169]);
print_mfilename(mfilename);

return;




function vs_iccpairs_plot_info_diff_pairs(infopairs)


edges = linspace(0,10,21);
xlabel1 = {'Information Difference (bits)'};
xlabel2 = {'Information Difference (bits)'};

ymaxdiff = [0.41 0.65];
ytick1 = 0:0.1:1; 
xtick2 = {[0.1 1 10 100]}; 
ytick2 = {[0.01 0.1 1]}; 

[larger,smaller] = vs_largersmaller(infopairs(:,1), infopairs(:,2));
datadiff = abs( larger - smaller );
data_median = nanmedian(datadiff);
data_mean = nanmean(datadiff);


subplot(1,2,1);
n = histc(datadiff, edges);
n = n ./ sum(n);
hb = bar(edges,n,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
hold on;
box off;
tickpref;
xtick = edges(1:2:end);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick1, 'yticklabel', ytick1);
range = max(edges) - min(edges);
xlim([min(edges)-0.025*range max(edges)+0.025*range]);
ylim([0 1.05*max(n)]);
ylabel('Proportion');
xlabel('Information Difference (bits)');


% Now calculate what we would expect if the values were randomly assigned.
nreps = 1000;
rand_median = zeros(1,nreps);
rand_mean = zeros(1,nreps);

for i = 1:nreps
   index1 = randperm(length(larger));
   index2 = randperm(length(larger));
   temp1 = larger(index1);
   temp2 = smaller(index2);
   datadiff = abs(temp1 - temp2);
   rand_median(i) = nanmedian(datadiff);
   rand_mean(i) = nanmean(datadiff);
end % (for i)

index = find(rand_median < data_median);
pval = length(index) / nreps;


% f1 = figure;
% hist(rand_median,25);
% min(rand_median)
% max(rand_median)
% pause
% close(f1);

subplot(1,2,2);
hold on;
[n,x] = hist(rand_median,50);
edges = center2edge(x);
n = histc(rand_median, edges);
n = n ./ sum(n);
hb = bar(edges,n,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
set(hb, 'edgecolor', 'none');
tickpref;
minmin = min([data_median min(rand_median)]);
maxmax = max([data_median max(rand_median)]);
range = maxmax - minmin;
xlim([minmin-0.05*range maxmax+0.05*range]);
ylim([0 max(n)]);
plot([data_median data_median], ylim, 'k-', 'linewidth', 2);
xlabel('Median Info Difference (bits)');
title(sprintf('pval = %.4f, n = %.0f', pval, nreps));
set(gcf,'position', [577 171 438 169]);
print_mfilename(mfilename);

return;








function vs_iccpairs_plot_rbar_pairs(rbarpairs)

[larger, smaller] = vs_largersmaller(rbarpairs(:,1), rbarpairs(:,2));
% scatter(larger, smaller, 25, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold on;
plot(larger, smaller, 'ko', 'markersize', 5);
tickpref;
box off;
xlabel('Ave Firing Rate Neuron 1 (sp/s)');
ylabel('Ave Firing Rate Neuron 2 (sp/s)');
set(gca,'xscale', 'log', 'yscale', 'log');

ytick = [0.1 1 10 50];
ylim([min(ytick) max(ytick)]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');

xtick = ytick;
xlim([min(xtick) max(xtick)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'xscale', 'log');

plot(xlim, ylim, 'k-');
text(min(xlim), 1.5*max(ylim), 'A');

index = ~isnan(larger) & ~isnan(smaller);
larger = larger(index);
smaller = smaller(index);
datadiff = abs(larger - smaller);
fprintf('\n');
fprintf('Rbar pairs\n');
[r,p] = corrcoef(log10(larger), log10(smaller));
fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
fprintf('Rbar pairs difference\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', length(datadiff),mean(datadiff), std(datadiff));
fprintf('MD = %.4f, MAD = %.4f, MX = %.4f\n', median(datadiff), mad(datadiff), max(datadiff));
fprintf('\n\n');


return;



function vs_iccpairs_plot_dt_pairs(dtpairs)
dtpairs = 1000 * dtpairs;

[larger, smaller] = vs_largersmaller(dtpairs(:,1), dtpairs(:,2));
% scatter(larger, smaller, 25, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold on;
plot(larger, smaller, 'ko', 'markersize', 5);
tickpref;
box off;
xlabel('Bin Size Neuron 1 (ms)');
ylabel('Bin Size Neuron 2 (ms)');

ytick = [0:5:25];
ylim([min(ytick) max(ytick)]);
set(gca,'ytick', ytick, 'yticklabel', ytick);

xtick = ytick;
xlim([min(xtick) max(xtick)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);

plot(xlim, ylim, 'k-');
text(min(xlim), 1.5*max(ylim), 'B');

index = ~isnan(larger) & ~isnan(smaller);
larger = larger(index);
smaller = smaller(index);
datadiff = abs(larger - smaller);
fprintf('\n');
fprintf('DT pairs\n');
[r,p] = corrcoef(larger, smaller);
fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
fprintf('Rbar pairs difference\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', length(datadiff),mean(datadiff), std(datadiff));
fprintf('MD = %.4f, MAD = %.4f, MX = %.4f\n', median(datadiff), mad(datadiff), max(datadiff));
fprintf('\n\n');


return;



function vs_iccpairs_plot_info_pairs(infopairs)

[larger, smaller] = vs_largersmaller(infopairs(:,1), infopairs(:,2));

% scatter(larger, smaller, 25, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold on;
plot(larger, smaller, 'ko', 'markersize', 5);
tickpref;
box off;
xlabel('Information Neuron 1 (bits/sp)');
ylabel('Information Neuron 2 (bits/sp)');

ytick = [0.2 1 10];
ylim([min(ytick) max(ytick)]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');

xtick = [0.2 1 10];
xlim([min(xtick) max(xtick)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'xscale', 'log');

plot(xlim, ylim, 'k-');
text(min(xlim), 1.5*max(ylim), 'C');

index = ~isnan(larger) & ~isnan(smaller);
larger = larger(index);
smaller = smaller(index);
datadiff = abs(larger - smaller);
fprintf('\n');
fprintf('Info pairs\n');
[r,p] = corrcoef(log10(larger), log10(smaller));
fprintf('r=%.4f,p=%.4f\n', r(2),p(2));
fprintf('Info pairs difference\n');
fprintf('N = %.0f, MN = %.4f, SD = %.4f\n', length(datadiff),mean(datadiff), std(datadiff));
fprintf('MD = %.4f, MAD = %.4f, MX = %.4f\n', median(datadiff), mad(datadiff), max(datadiff));
fprintf('\n\n');


return;




















