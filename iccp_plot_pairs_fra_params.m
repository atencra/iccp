function [fr, cf, bw, bwlog, Q, latency] = iccp_plot_pairs_fra_params(fr, cf, bw, bwlog, Q, latency)
% iccp_plot_pairs_fra_params Pairwise/difference plots for FRA parameters
% 
%     iccp_plot_pairs_fra_params; looks for '*-fracmb-pairs-params.mat' in
%     the current directory, extracts the data, and compares FRA values for 
%     pairs of neurons. Scatter plots for the values, and difference histograms 
%     for each pair, are shown. The parameters that are plotted are firing rate,
%     CF, BW, log BW, Q, and latency.
% 
%     [fr, cf, bw, bwlog, Q, latency] = iccp_plot_pairs_fra_params; returns the
%     data used to make the plots.
% 
% 
%     iccp_plot_pairs_fra_params(fr, cf, bw, bwlog, Q, latency); uses the 
%     previously returned data to quickly make the plots. This is helpful
%     because it doesn't read through all the files in the directory.
%
%     You can indicate specific plots by using the following:
%
%     iccp_plot_pairs_fra_params(str); where str is a string specifying
%     the type of plot. str may be 'fr', 'cf', 'bw', 'q', or 'latency'.
%     Data is read from the files in this form.
%
%     iccp_plot_pairs_fra_params(str, data); make a plot without searching
%     through the files. The str takes on values as before, and data is an
%     Nx2 array that corresponds to str.



narginchk(0,6);

if ( nargin == 0 )
   [position, fr, cf, bw, bwlog, Q, latency] = get_pairs_fra_params;
end

close all;

% vs_plot_pairs_fra_fr_scatter_hist(fr);

% vs_plot_pairs_fra_cf_scatter_hist(cf);



% vs_plot_pairs_fra_bw_scatter_hist(bw);

% vs_plot_pairs_fra_bw_diff_rand(bw);



vs_plot_pairs_fra_q_scatter_hist(Q);

% vs_plot_pairs_fra_q_diff_rand(Q)


% vs_plot_pairs_fra_latency_scatter_hist(latency);
% 
% plot_pairs_population_diff_layers_fra_params(fr, cf, bw, bwlog, Q, latency);


% label = {'BW10', 'BW20', 'BW30', 'BW40'};
% for i = 1:4
%    data = [bwlog(:,i) bwlog(:,i+4)];
%    index = ~isnan(data(:,1)) & ~isnan(data(:,2));
%    data = data(index,:);
%    r = corrcoef(data(:,1), data(:,2));
%    r = r(2);
%    bwdata{i} = data;
%    bwr{i} = r;
%    bwdiff{i} = abs(data(:,1) - data(:,2));
%    bwdiff_oct{i} = abs( log2(data(:,1)./data(:,2)) );
% 
%    datadiff = bwdiff{i};
%    fprintf('\n%s\n',label{i});
% %    [r,p] = corrcoef(larger, smaller);
% %    fprintf('r = %.4f, p = %.4f\n', r(2), p(2) );
%    fprintf('MN = %.4f, MD = %.4f, SD = %.4f, SE = %.4f, n = %.0f\n',...
%       mean(datadiff), median(datadiff), std(datadiff), ...
%       std(datadiff)./sqrt( length(datadiff) ), length(datadiff) );
%    fprintf('\n');
% 
% end % (for i)


return;



function [position, fr, cf, bw, bwlog, Q, latency] = get_pairs_fra_params

position = [];
fr = [];
cf = [];
bw1 = [];
bw2 = [];
BWL1 = [];
BWL2  = [];
Q = [];
latency = [];

% d = dir('*-fra-params.mat');
d = dir('*-fracmb-pairs-params.mat');


for n = 1:length(d)
    
    filename = d(n).name;
    s = load(filename, 'params');
    params = s.params;
    
    chan = [params.chan];
    chan_unique = unique(chan);
    
    for i = 1:length(chan_unique)
        
        index = find( chan_unique(i) == chan );
        
        if length(index)>1
            
            cmb = nchoosek(index, 2); % determine all possible pairwise combinations
            
            [nr, nc] = size(cmb);
            
            for j = 1:nr
                
                index1 = cmb(j,1);
                index2 = cmb(j,2);
                
                position = [position; params(index1).position];
                
                cf = [cf; params(index1).cf params(index2).cf];
                cf_1 = params(index1).cf;
                cf_2 = params(index2).cf;
                
                %Neuron 1 Bandwidths
                if isempty(params(index1).bw10)
                    bw10_1 = NaN;
                    BWL10_1 = NaN;
                else
                    bw10_1 = params(index1).bw10(2)-params(index1).bw10(1);
                    BWL10_1 = abs(log2(params(index1).bw10(2))- log2(params(index1).bw10(1)));
                end
                
                if isempty(params(index1).bw20)
                    bw20_1 = NaN;
                    BWL20_1 = NaN;
                else
                    bw20_1 = params(index1).bw20(2)-params(index1).bw20(1);
                    BWL20_1 = abs(log2(params(index1).bw20(2))- log2(params(index1).bw20(1)));
                end
                
                if isempty(params(index1).bw30)
                    bw30_1 = NaN;
                    BWL30_1 = NaN;
                else
                    bw30_1 = params(index1).bw30(2)-params(index1).bw30(1);
                    BWL30_1 = abs(log2(params(index1).bw30(2))- log2(params(index1).bw30(1)));
                end
                
                if isempty(params(index1).bw40)
                    bw40_1 = NaN;
                    BWL40_1 = NaN;
                else
                    bw40_1 = params(index1).bw40(2)-params(index1).bw40(1);
                    BWL40_1 = abs(log2(params(index1).bw40(2))- log2(params(index1).bw40(1)));
                end
                bw1 = [bw1; bw10_1 bw20_1 bw30_1 bw40_1];
                BWL1 = [BWL1; BWL10_1 BWL20_1 BWL30_1 BWL40_1];


               
                %Neuron 2 Bandwidths
                if isempty(params(index2).bw10)
                    bw10_2 = NaN;
                    BWL10_2 = NaN;
                else
                    bw10_2 = params(index2).bw10(2)-params(index2).bw10(1);
                    BWL10_2 = abs(log2(params(index2).bw10(2))- log2(params(index2).bw10(1)));
                end
                
                if isempty(params(index2).bw20)
                    bw20_2 = NaN;
                    BWL20_2 = NaN;
                else
                    bw20_2 = params(index2).bw20(2)-params(index2).bw20(1);
                    BWL20_2 = abs(log2(params(index2).bw20(2))- log2(params(index2).bw20(1)));
                end
                
                if isempty(params(index2).bw30)
                    bw30_2 = NaN;
                    BWL30_2 = NaN;
                else
                    bw30_2 = params(index2).bw30(2)-params(index2).bw30(1);
                    BWL30_2 = abs(log2(params(index2).bw30(2))- log2(params(index2).bw30(1)));
                end
                
                if isempty(params(index2).bw40)
                    bw40_2 = NaN;
                    BWL40_2 = NaN;
                else
                    bw40_2 = params(index2).bw40(2)-params(index2).bw40(1);
                    BWL40_2 = abs(log2(params(index2).bw40(2))- log2(params(index2).bw40(1)));
                end
                bw2 = [bw2; bw10_2 bw20_2 bw30_2 bw40_2];
                BWL2 = [BWL2; BWL10_2 BWL20_2 BWL30_2 BWL40_2];
                
%               Q Values
                Q10 = [cf_1./bw10_1 cf_2./bw10_2];
                Q20 = [cf_1./bw20_1 cf_2./bw20_2];
                Q30 = [cf_1./bw30_1 cf_2./bw30_2];
                Q40 = [cf_1./bw40_1 cf_2./bw40_2];
                Q = [Q; Q10 Q20 Q30 Q40];                
                
                % Latency and Firing Rate
                latency = [latency; params(index1).latency params(index2).latency];
                
                fr = [fr; max(max(params(index1).rlcurve(:,2:6))) max(max(params(index2).rlcurve(:,2:6)))];
                
            end % (for j)
            
        end % (if/else)
        
    end % (for i)
    
    bw = [bw1 bw2];
    bwlog = [BWL1 BWL2];
    
end % (for n)
return;




function vs_plot_pairs_fra_cf_scatter_hist(cf)

% CF
cfdiff = abs( log2( cf(:,2)./cf(:,1) ) );
dstats = simple_stats(cfdiff);

cftick = [0.125 0.25 0.5 1 2 4 8 16 32];
figure;
subplot(1,2,1);
hold on
axis([0.2 32 0.2 32])
plot(xlim, ylim, 'k-');


[larger, smaller] = iccp_largersmaller(cf(:,1), cf(:,2));

[absc, ord] = iccp_randomize_columns(cf(:,1), cf(:,2));


hold on;

%scatter(larger,smaller, 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6]);

scatter(absc, ord, 10, 'MarkerEdgeColor', 'black'); %, 'MarkerFaceColor', [0.8 0.8 0.8]);


tickpref;
box off;
xlabel('CF (kHz) Neuron 1')
ylabel('CF (kHz) Neuron 2')
set(gca,'xtick', cftick, 'xticklabel', cftick);
set(gca,'ytick', cftick, 'yticklabel', cftick);
set(gca,'xscale', 'log', 'yscale', 'log');

subplot(1,2,2);
edges = linspace(0,3,16);
n = histc(cfdiff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
xlabel('CF Difference (oct)')
ylabel('Proportion');
tickpref
box off
set(gca,'ytick', 0:0.25:0.75, 'yticklabel', 0:0.25:0.75);
set(gca,'xtick', 0:1:4, 'xticklabel', 0:1:4);
ylim([0 0.75])
xlim([min(edges) max(edges)]);

set(gcf,'position', [1255 657 438 161]);
print_mfilename(mfilename);

[r,p] = corrcoef(log10(larger), log10(smaller) );
fprintf('\n');
fprintf('CF Comparison\n');
fprintf('r = %.4f, p = %.4f\n',r(2), p(2) );

fprintf('\n');
fprintf('CF Difference (oct)\n');
fprintf('MN = %.4f, SD = %.4f, SE = %.4f, MD = %.4f, MAD = %.4f, n = %.0f\n',...
   dstats.mn, dstats.sd, dstats.se, dstats.md, dstats.mad, length(cf(:,1)) );



logtransform = 1;
[rpop_med, rpop_ci, pval] = ...
    iccp_pairwise_corr_rand_test(cf(:,1), cf(:,2), logtransform);

fprintf('\n');
fprintf('CF pairs - randomized\n');
fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));

fprintf('\n\n');




return;




function vs_plot_pairs_fra_bw_scatter_hist(bw)



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
ytick = {0:0.05:0.25, 0:0.05:0.35, 0:0.05:0.35, 0:0.05:0.35}; 
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

    [larger, smaller] = iccp_largersmaller(bwdata{i}(:,1),bwdata{i}(:,2));

    [absc, ord] = iccp_randomize_columns(bwdata{i}(:,1),bwdata{i}(:,2));

%     plot(larger,smaller, 'ko', ...
%         'MarkerEdgeColor', 'black', ...
%         'MarkerSize', 2);


    plot(absc, ord, 'ko', ...
        'MarkerEdgeColor', 'black', ...
        'MarkerSize', 2);


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



    logtransform = 0;
    [rpop_med, rpop_ci, pval] = ...
        iccp_pairwise_corr_rand_test(bwdata{i}(:,1), bwdata{i}(:,2), logtransform);
    fprintf('\n');
    fprintf(sprintf('%s pairs - randomized\n',label{i}));
    fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
    fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
    fprintf('\n\n');

end % (for i)


set(gcf,'Position',[400 100 321 551])
print_mfilename(mfilename);


return;




function vs_plot_pairs_fra_bw_diff_rand(bw)

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
ytick = {0:0.05:0.25, 0:0.05:0.35, 0:0.05:0.35, 0:0.05:0.35}; 
minmin = 0.125;
maxmax = 32;

for i = 1:4

   data = bwdata{i};
   datadiff = bwdiff_oct{i};
   data_median = nanmedian(datadiff);
   data_mean = nanmean(datadiff);

   [larger, smaller] = vs_largersmaller(bwdata{i}(:,1),bwdata{i}(:,2));

   subplot(4,2,(i-1)*2+1);
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
   fprintf('MN = %.4f, MD = %.4f, SD = %.4f, SE = %.4f, n = %.0f\n',...
      mean(datadiff), median(datadiff), std(datadiff), ...
      std(datadiff)./sqrt( length(datadiff) ), length(datadiff) );
   fprintf('\n');






   % Now calculate what we would expect if the values were randomly assigned.
   nreps = 1000;
   rand_median = zeros(1,nreps);
   rand_mean = zeros(1,nreps);

   for ii = 1:nreps
      index1 = randperm(length(larger));
      index2 = randperm(length(larger));
      temp1 = larger(index1);
      temp2 = smaller(index2);
      randdiff = abs( log2( temp1 ./ temp2 ) );
      rand_median(ii) = nanmedian(randdiff);
      rand_mean(ii) = nanmean(randdiff);
   end % (for i)

   index = find(rand_median < data_median);
   pval = length(index) / nreps;



   subplot(4,2,(i)*2);
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
   xlabel('Median Bin Size Diff');
   title(sprintf('pval = %.4f, n = %.0f', pval, nreps));

end % (for i)


set(gcf,'Position',[400 100 321 551])
print_mfilename(mfilename);


return;




function vs_plot_pairs_fra_q_scatter_hist(Q)


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

ydiffmax = [0.25 0.3 0.35 0.35];
ytick = {0:0.05:0.25, 0:0.05:0.3, 0:0.05:0.35, 0:0.05:0.35}; 
markersize = 1;

for i = 1:4

   data = qdata{i};
   datadiff = qdiff_oct{i};

   subplot(4,2,(i-1)*2+1);
   hold on
   xlim([minmin maxmax]);
   ylim([minmin maxmax]);
   plot(xlim, ylim, 'k-');

%    [larger, smaller] = iccp_largersmaller(qdata{i}(:,1),qdata{i}(:,2));
% 
%     plot(larger,smaller, 'ko', ...
%         'MarkerEdgeColor', 'black', ...
%         'MarkerSize', 2);


    [absc, ord] = iccp_randomize_columns(qdata{i}(:,1),qdata{i}(:,2));

    plot(absc, ord, 'ko', ...
        'MarkerEdgeColor', 'black', ...
        'MarkerSize', 2);




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


%    fprintf('\n%s\n',label{i});
%    [r,p] = corrcoef(larger, smaller);
%    fprintf('r = %.4f, p = %.4f\n', r(2), p(2) );
% 
%    dstats = simple_stats(datadiff);
% 
%    fprintf('\n');
%    fprintf('MN = %.4f, SD = %.4f, SE = %.4f, MD = %.4f, MAD = %.4f, n = %.0f\n',...
%       dstats.mn, dstats.sd, dstats.se, dstats.md, dstats.mad, length(datadiff) );
%    fprintf('\n');


    logtransform = 0;
    [rpop_med, rpop_ci, pval] = ...
        iccp_pairwise_corr_rand_test(qdata{i}(:,1), qdata{i}(:,2), logtransform);
    fprintf('\n');
    fprintf(sprintf('%s pairs - randomized\n',label{i}));
    fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
    fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
    fprintf('\n\n');

end % (for i)

set(gcf,'Position',[400 100 321 551])
print_mfilename(mfilename);

return;




function vs_plot_pairs_fra_q_diff_rand(Q)


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

ydiffmax = [0.25 0.3 0.35 0.35];
ytick = {0:0.05:0.25, 0:0.05:0.3, 0:0.05:0.35, 0:0.05:0.35}; 
markersize = 1;

for i = 1:4

   data = qdata{i};
   datadiff = qdiff_oct{i};
   data_median = nanmedian(datadiff);
   data_mean = nanmean(datadiff);
   [larger, smaller] = vs_largersmaller(qdata{i}(:,1),qdata{i}(:,2));

   subplot(4,2,(i-1)*2+1);
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
   fprintf('MN = %.4f, MD = %.4f, SD = %.4f, SE = %.4f, n = %.0f\n',...
      mean(datadiff), median(datadiff), std(datadiff), ...
      std(datadiff)./sqrt( length(datadiff) ), length(datadiff) );
   fprintf('\n');


   % Now calculate what we would expect if the values were randomly assigned.
   nreps = 1000;
   rand_median = zeros(1,nreps);
   rand_mean = zeros(1,nreps);

   for ii = 1:nreps
      index1 = randperm(length(larger));
      index2 = randperm(length(larger));
      temp1 = larger(index1);
      temp2 = smaller(index2);
      randdiff = abs( log2( temp1 ./ temp2 ) );
      rand_median(ii) = nanmedian(randdiff);
      rand_mean(ii) = nanmean(randdiff);
   end % (for i)

   index = find(rand_median < data_median);
   pval = length(index) / nreps;


   subplot(4,2,(i)*2);
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
   xlabel('Median Bin Size Diff');
   title(sprintf('pval = %.4f, n = %.0f', pval, nreps));

end % (for i)

set(gcf,'Position',[400 100 321 551])
print_mfilename(mfilename);

return;




function vs_plot_pairs_fra_latency_scatter_hist(latency)

% Latency 

lat_diff = abs(latency(:,2) - latency(:,1));
med_lat_diff = median(lat_diff);
mean_lat_diff = mean(lat_diff);
std_lat_diff = std(lat_diff);
R_lat = corrcoef(latency(:,1),latency(:,2));
r_lat = R_lat(2);


figure;
subplot(2,1,1);
hold on
axis([0 40 0 40])
plot(xlim, ylim, 'k-');
[larger, smaller] = vs_largersmaller(latency(:,1), latency(:,2));
scatter(larger,smaller,25, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
tickpref
box off
xlabel('Latency of Neuron 1')
ylabel('Latency of Neuron 2')
title(sprintf('r = %.2f',r_lat))
axis([0 40 0 40])
suptitle('Latency')

subplot(2,1,2);
edges = linspace(0,10,21);
n = histc(lat_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('MN= %.2f, MD= %.2f, SD= %.2f',mean_lat_diff, med_lat_diff,std_lat_diff))
xlabel('Latency Difference (ms)')
tickpref
box off
set(gca,'ytick', 0:0.1:0.3, 'yticklabel', 0:0.1:0.3);
ylim([0 0.3])
xlim([min(edges) max(edges)]);
set(gcf,'Position',[177   485   321   404]);
print_mfilename(mfilename);

return;



function plot_pairs_population_diff_layers_fra_params(fr, cf, bw, bwlog, Q, latency)

% Firing Rate
fr_diff = abs(fr(:,2) - fr(:,1));
med_fr_diff = median(fr_diff);
mean_fr_diff = mean(fr_diff);
std_fr_diff = std(fr_diff);
R_fr = corrcoef(fr(:,1),fr(:,2));
r_fr = R_fr(2);

figure
set(gcf,'Position',[100 100 700 700])
subplot(5,2,[1 6]);
hold on
plot([0 100],[0 100],'k','LineWidth',1.5)
scatter(fr(:,1),fr(:,2),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold off
tickpref
box off
xlabel('Firing Rate of Neuron 1')
ylabel('Firing Rate of Neuron 2')
axis([0 100 0 100])
suptitle('Firing Rates')
title(sprintf('r = %.2f',r_fr))

subplot(5,2,[9 10]);
edges = linspace(0,60,20);
n = histc(fr_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Firing Rate Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_fr_diff, med_fr_diff,std_fr_diff))
xlabel('Spikes Per Second')
tickpref
box off
ylim([0 1.05*max(n)])
print_mfilename(mfilename);

% CF
CF_log_diff = abs(log2(cf(:,2) - log2(cf(:,1))));
med_CF_diff = median(CF_log_diff);
mean_CF_diff = mean(CF_log_diff);
std_CF_diff = std(CF_log_diff);
R_CF = corrcoef(cf(:,1),cf(:,2));
r_CF = R_CF(2);












figure
set(gcf,'Position',[100 100 700 700])
subplot(5,2,[1 6]);
hold on
plot([min(cf(:,1))*0.95 max(cf(:,1))*1.05],[min(cf(:,2))*0.95 max(cf(:,2))*1.05],'k','LineWidth', 1.5)
scatter(cf(:,1),cf(:,2),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold off
tickpref
box off
xlabel('CF of Neuron 1')
ylabel('CF of Neuron 2')
axis([0 30 0 30])
set(gca,'xscale', 'log', 'yscale', 'log');
suptitle('Characteristic Frequency')
title(sprintf('CF: r = %.2f',r_CF))

subplot(5,2,[9 10]);
edges = linspace(0,4,21);
n = histc(CF_log_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('CF Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_CF_diff, med_CF_diff,std_CF_diff))
xlabel('Octaves')
tickpref
box off
ylim([0 1.05*max(n)])
print_mfilename(mfilename);


% BW10 
% figure
% subplot(5,2,[1 6]);
% scatter(bw(:,1),bw(:,5),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW10 of Neuron 1')
% ylabel('BW10 of Neuron 2')
% axis([0 15 0 15])
% suptitle('Bandwidth at 10 dB SPL above Threshold')

% Find and remove NaN entries to calculate R-value 
BW10 = [bw(:,1) bw(:,5)];
L = length(BW10);
locations = [];
for i = 1:L
    if isnan(BW10(i,1)) || isnan(BW10(i,2))
        locations = [locations i];
    end
end
BW10(locations,:) = [];
R_bw10 = corrcoef(BW10(:,1),BW10(:,2));
r_bw10 = R_bw10(2);

BW10_log_diff = abs(log2(BW10(:,2) - log2(BW10(:,1))));
med_BW10_diff = median(BW10_log_diff);
mean_BW10_diff = mean(BW10_log_diff);
std_BW10_diff = std(BW10_log_diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,3,21);
% n = histc(BW10_log_diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('BW10 Difference between Neuron 1 and Neuron 2')
% tickpref
% box off

% figure
% scatter(bw(:,1),bw(:,5),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% set(gca,'xscale', 'log', 'yscale', 'log');
% suptitle('Bandwidth at 10 dB SPL above Threshold')
% ylabel('BW10 of Neuron 2')
% xlabel('BW10 of Neuron 1')

% print_mfilename('plot_pairs_fra_params');

% BW20
% figure
% subplot(5,2,[1 6]);
% scatter(bw(:,2),bw(:,6),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW20 of Neuron 1')
% ylabel('BW20 of Neuron 2')
% axis([0 27 0 27])
% suptitle('Bandwidth at 20 dB SPL above Threshold')

% Find and remove NaN entries
BW20 = [bw(:,2) bw(:,6)];
L = length(BW20);
locations = [];
for i = 1:L
    if isnan(BW20(i,1)) || isnan(BW20(i,2))
        locations = [locations i];
    end
end
BW20(locations,:) = [];
R_bw20 = corrcoef(BW20(:,1),BW20(:,2));
r_bw20 = R_bw20(2);

BW20_log_diff = abs(log2(BW20(:,2) - log2(BW20(:,1))));
med_BW20_diff = median(BW20_log_diff);
mean_BW20_diff = mean(BW20_log_diff);
std_BW20_diff = std(BW20_log_diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,3,21);
% n = histc(BW20_log_diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('BW20 Difference between Neuron 1 and Neuron 2')
% tickpref
% box off

% figure
% scatter(bw(:,2),bw(:,6),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% set(gca,'xscale', 'log', 'yscale', 'log');
% title('Bandwidth at 20 dB SPL above Threshold')
% ylabel('BW10 of Neuron 2')
% xlabel('BW10 of Neuron 1')

% print_mfilename('plot_pairs_fra_params');

% BW30
% figure
% subplot(5,2,[1 6]);
% scatter(bw(:,3),bw(:,7),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW30 of Neuron 1')
% ylabel('BW30 of Neuron 2')
% axis([0 35 0 35])
% suptitle('Bandwidth at 30 dB SPL above Threshold')

% Find and remove NaN entries to calculate R-value 
BW30 = [bw(:,3) bw(:,7)];
L = length(BW30);
locations = [];
for i = 1:L
    if isnan(BW30(i,1)) || isnan(BW30(i,2))
        locations = [locations i];
    end
end
BW30(locations,:) = [];
R_bw30 = corrcoef(BW30(:,1),BW30(:,2));
r_bw30 = R_bw30(2);

BW30_log_diff = abs(log2(BW30(:,2) - log2(BW30(:,1))));
med_BW30_diff = median(BW30_log_diff);
mean_BW30_diff = mean(BW30_log_diff);
std_BW30_diff = std(BW30_log_diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,3,21);
% n = histc(BW30_log_diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('BW30 Difference between Neuron 1 and Neuron 2')
% tickpref
% box off

% figure
% scatter(bw(:,3),bw(:,7),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% set(gca,'xscale', 'log', 'yscale', 'log');
% title('Bandwidth at 30 dB SPL above Threshold')
% ylabel('BW10 of Neuron 2')
% xlabel('BW10 of Neuron 1')

% print_mfilename('plot_pairs_fra_params');

% BW40
% figure
% subplot(5,2,[1 6]);
% scatter(bw(:,4),bw(:,8),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW40 of Neuron 1')
% ylabel('BW40 of Neuron 2')
% axis([0 40 0 40])
% suptitle('Bandwidth at 40 dB SPL above Threshold')

% Find and remove NaN entries to calculate R-value 
BW40 = [bw(:,4) bw(:,8)];
L = length(BW40);
locations = [];
for i = 1:L
    if isnan(BW40(i,1)) || isnan(BW40(i,2))
        locations = [locations i];
    end
end
BW40(locations,:) = [];
R_bw40 = corrcoef(BW40(:,1),BW40(:,2));
r_bw40 = R_bw40(2);

BW40_log_diff = abs(log2(BW40(:,2) - log2(BW40(:,1))));
med_BW40_diff = median(BW40_log_diff);
mean_BW40_diff = mean(BW40_log_diff);
std_BW40_diff = std(BW40_log_diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,4,21);
% n = histc(BW40_log_diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('BW40 Difference between Neuron 1 and Neuron 2')
% tickpref
% box off

% figure
% scatter(bw(:,4),bw(:,8),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% set(gca,'xscale', 'log', 'yscale', 'log');
% title('Bandwidth at 40 dB SPL above Threshold')
% ylabel('BW10 of Neuron 2')
% xlabel('BW10 of Neuron 1')

% print_mfilename('plot_pairs_fra_params');

% Plot of BW
figure
set(gcf,'Position',[100 100 1400 700])
suptitle('Bandwidth')

% subplot(7,8,[1 2]);
% scatter(bw(:,1),bw(:,5),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW10 of Neuron 1')
% ylabel('BW10 of Neuron 2')
% axis([0 15 0 15])
% title('BW10')

subplot(7,4,[1 2]);
hold on
plot([min(bw(:,1))*0.95 max(bw(:,1))*1.05],[min(bw(:,5))*0.95 max(bw(:,5))*1.05],'k','LineWidth', 1.5)
scatter(bw(:,1),bw(:,5),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([min(bw(:,1))*0.95 max(bw(:,1))*1.05])
ylim([min(bw(:,5))*0.95 max(bw(:,5))*1.05])
hold off
tickpref
box off
set(gca,'xscale', 'log', 'yscale', 'log');
ylabel('BW10 of Neuron 2')
xlabel('BW10 of Neuron 1')
title(sprintf('BW10: r = %.2f',r_bw10))

subplot(7,4,[3 4]);
edges = linspace(0,3,21);
n = histc(BW10_log_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('BW10 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_BW10_diff, med_BW10_diff,std_BW10_diff))
xlabel('Octaves')
tickpref
box off

% subplot(7,8,[17 18]);
% scatter(bw(:,2),bw(:,6),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW20 of Neuron 1')
% ylabel('BW20 of Neuron 2')
% axis([0 27 0 27])
% title('Bandwidth at 20 dB SPL above Threshold')

subplot(7,4,[9 10]);
hold on
plot([min(bw(:,2))*0.95 max(bw(:,2))*1.05],[min(bw(:,6))*0.95 max(bw(:,6))*1.05],'k','LineWidth', 1.5)
scatter(bw(:,2),bw(:,6),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([min(bw(:,2))*0.95 max(bw(:,6))*1.05])
ylim([min(bw(:,2))*0.95 max(bw(:,6))*1.05])
tickpref
box off
set(gca,'xscale', 'log', 'yscale', 'log');
title(sprintf('BW20: r = %.2f',r_bw20))
ylabel('BW20 of Neuron 2')
xlabel('BW20 of Neuron 1')

subplot(7,4,[11 12]);
edges = linspace(0,3,21);
n = histc(BW20_log_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('BW20 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_BW20_diff, med_BW20_diff,std_BW20_diff))
xlabel('Octaves')
tickpref
box off

% subplot(7,4,[17 18]);
% scatter(bw(:,3),bw(:,7),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW30 of Neuron 1')
% ylabel('BW30 of Neuron 2')
% axis([0 35 0 35])
% title('Bandwidth at 30 dB SPL above Threshold')

subplot(7,4,[17 18]);
hold on
plot([min(bw(:,3))*0.95 max(bw(:,3))*1.05],[min(bw(:,7))*0.95 max(bw(:,7))*1.05],'k','LineWidth', 1.5)
scatter(bw(:,3),bw(:,7),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([min(bw(:,3))*0.95 max(bw(:,7))*1.05])
ylim([min(bw(:,7))*0.95 max(bw(:,7))*1.05])
hold off
tickpref
box off
set(gca,'xscale', 'log', 'yscale', 'log');
title(sprintf('BW30: r = %.2f',r_bw30))
ylabel('BW30 of Neuron 2')
xlabel('BW30 of Neuron 1')

subplot(7,4,[19 20]);
edges = linspace(0,3,21);
n = histc(BW30_log_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('BW30 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_BW30_diff, med_BW30_diff,std_BW30_diff))
xlabel('Octaves')
tickpref
box off

% subplot(7,8,[49 50]);
% scatter(bw(:,4),bw(:,8),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% box off
% xlabel('BW40 of Neuron 1')
% ylabel('BW40 of Neuron 2')
% axis([0 40 0 40])
% title('Bandwidth at 40 dB SPL above Threshold')

subplot(7,4,[25 26]);
hold on
plot([min(bw(:,4))*0.95 max(bw(:,4))*1.05],[min(bw(:,8))*0.95 max(bw(:,8))*1.05],'k','LineWidth', 1.5)
scatter(bw(:,4),bw(:,8),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([min(bw(:,4))*0.95 max(bw(:,4))*1.05])
ylim([min(bw(:,8))*0.95 max(bw(:,8))*1.05])
hold off
tickpref
box off
set(gca,'xscale', 'log', 'yscale', 'log');
title(sprintf('BW40: r = %.2f',r_bw40))
ylabel('BW40 of Neuron 2')
xlabel('BW40 of Neuron 1')

subplot(7,4,[27 28]);
edges = linspace(0,4,21);
n = histc(BW40_log_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('BW40 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_BW40_diff, med_BW40_diff,std_BW40_diff))
xlabel('Octaves')
tickpref
box off

print_mfilename(mfilename);



















% Q10
% figure
% subplot(5,2,[1 6]);
% scatter(Q(:,1),Q(:,5),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% xlabel('Q10 of Neuron 1')
% ylabel('Q10 of Neuron 2')
% axis([0 8 0 8])
% set(gca,'xscale', 'log', 'yscale', 'log');
% suptitle('Q at 10 dB SPL above Threshold')

% Find and remove NaN entries to calculate R-value 
Q10 = [Q(:,1) Q(:,2)];
L = length(Q10);
locations = [];
for i = 1:L
    if isnan(Q10(i,1)) || isnan(Q10(i,2))
        locations = [locations i];
    end
end
Q10(locations,:) = [];
R_Q10 = corrcoef(Q10(:,1),Q10(:,2));
r_Q10 = R_Q10(2);

Q_10diff = abs(Q10(:,2) - Q10(:,1)); 
med_Q10_diff = median(Q_10diff);
mean_Q10_diff = mean(Q_10diff);
std_Q10_diff = std(Q_10diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,4,21);
% n = histc(Q_10diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('Q10 Difference between Neuron 1 and Neuron 2, n=724')
% tickpref
% box off
% ylim([0 1.05*max(n)])

% print_mfilename('plot_pairs_fra_params');

% Q20
% figure
% subplot(5,2,[1 6]);
% scatter(Q(:,2),Q(:,6),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% xlabel('Q20 of Neuron 1')
% ylabel('Q20 of Neuron 2')
% axis([0 8 0 8])
% set(gca,'xscale', 'log', 'yscale', 'log');
% suptitle('Q at 20 dB SPL above Threshold')

% Find and remove NaN entries to calculate R-value 
Q20 = [Q(:,3) Q(:,4)];
L = length(Q20);
locations = [];
for i = 1:L
    if isnan(Q20(i,1)) || isnan(Q20(i,2))
        locations = [locations i];
    end
end
Q20(locations,:) = [];
R_Q20 = corrcoef(Q20(:,1),Q20(:,2));
r_Q20 = R_Q20(2);

Q_20diff = abs(Q20(:,2) - Q20(:,1)); 
med_Q20_diff = median(Q_20diff);
mean_Q20_diff = mean(Q_20diff);
std_Q20_diff = std(Q_20diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,4,21);
% n = histc(Q_20diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('Q20 Difference between Neuron 1 and Neuron 2, n=724')
% tickpref
% box off
% ylim([0 1.05*max(n)])

% print_mfilename('plot_pairs_fra_params');

% Q30
% figure
% subplot(5,2,[1 6]);
% scatter(Q(:,3),Q(:,7),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% xlabel('Q30 of Neuron 1')
% ylabel('Q30 of Neuron 2')
% axis([0 6 0 6])
% set(gca,'xscale', 'log', 'yscale', 'log');
% suptitle('Q at 30 dB SPL above Threshold')

% Find and remove NaN entries to calculate R-value 
Q30 = [Q(:,5) Q(:,6)];
L = length(Q30);
locations = [];
for i = 1:L
    if isnan(Q30(i,1)) || isnan(Q30(i,2))
        locations = [locations i];
    end
end
Q30(locations,:) = [];
R_Q30 = corrcoef(Q30(:,1),Q30(:,2));
r_Q30 = R_Q30(2);

Q_30diff = abs(Q30(:,2) - Q30(:,1)); 
med_Q30_diff = median(Q_30diff);
mean_Q30_diff = mean(Q_30diff);
std_Q30_diff = std(Q_30diff);


% subplot(5,2,[9 10]);
% edges = linspace(0,1.5,21);
% n = histc(Q_30diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('Q30 Difference between Neuron 1 and Neuron 2, n=724')
% tickpref
% box off
% ylim([0 1.05*max(n)])

% print_mfilename('plot_pairs_fra_params');

% Q40
% figure
% subplot(5,2,[1 6]);
% scatter(Q(:,4),Q(:,8),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
% tickpref
% xlabel('Q40 of Neuron 1')
% ylabel('Q40 of Neuron 2')
% axis([0 4 0 4])
% set(gca,'xscale', 'log', 'yscale', 'log');
% suptitle('Q at 40 dB SPL above Threshold, R = ')

% Find and remove NaN entries to calculate R-value 
Q40 = [Q(:,7) Q(:,8)];
L = length(Q40);
locations = [];
for i = 1:L
    if isnan(Q40(i,1)) || isnan(Q40(i,2))
        locations = [locations i];
    end
end
Q40(locations,:) = [];
R_Q40 = corrcoef(Q40(:,1),Q40(:,2));
r_Q40 = R_Q40(2);

Q_40diff = abs(Q40(:,2) - Q40(:,1)); 
med_Q40_diff = median(Q_40diff);
mean_Q40_diff = mean(Q_40diff);
std_Q40_diff = std(Q_40diff);

% subplot(5,2,[9 10]);
% edges = linspace(0,1.5,21);
% n = histc(Q_40diff, edges);
% n = n./sum(n);
% hb = bar(edges,n,'histc');
% set(hb, 'FaceColor', [0.6 0.6 0.6])
% title('Q40 Difference between Neuron 1 and Neuron 2, n=724')
% tickpref
% box off
% ylim([0 1.05*max(n)])

% print_mfilename('plot_pairs_fra_params');

% Q Plots
figure 
set(gcf,'Position',[100 100 1400 700])
suptitle('Q')

subplot(7,4,[1 2]);
hold on
plot([0.05 10],[0.05 10], 'k')
scatter(Q(:,1),Q(:,2),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([0.05 10])
ylim([0.05 10])
hold off
tickpref
box off
xlabel('Q10 of Neuron 1')
ylabel('Q10 of Neuron 2')
title(sprintf('Q10: r = %.2f',r_Q10))
set(gca,'xscale', 'log', 'yscale', 'log');
xlim([min(Q(:,1))*0.95 max(Q(:,1))*1.05])
ylim([min(Q(:,2))*0.95 max(Q(:,2))*1.05])

subplot(7,4,[3 4]);
edges = linspace(0,4,21);
n = histc(Q_10diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Q10 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_Q10_diff, med_Q10_diff,std_Q10_diff))
xlabel('Octaves')
tickpref
box off
ylim([0 1.05*max(n)])

subplot(7,4,[9 10]);
hold on
plot([0.05 10],[0.05 10], 'k')
scatter(Q(:,3),Q(:,4),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([0.05 10])
ylim([0.05 10])
hold off
tickpref
box off
xlabel('Q20 of Neuron 1')
ylabel('Q20 of Neuron 2')
title(sprintf('Q20: r = %.2f',r_Q20))
set(gca,'xscale', 'log', 'yscale', 'log');
xlim([min(Q(:,3))*0.95 max(Q(:,3))*1.05])
ylim([min(Q(:,4))*0.95 max(Q(:,4))*1.05])

subplot(7,4,[11 12]);
edges = linspace(0,4,21);
n = histc(Q_20diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Q20 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_Q20_diff, med_Q20_diff,std_Q20_diff))
xlabel('Octaves')
tickpref
box off
ylim([0 1.05*max(n)])

subplot(7,4,[17 18]);
hold on
plot([0.05 10],[0.05 10], 'k')
scatter(Q(:,5),Q(:,6),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([0.05 10])
ylim([0.05 10])
hold off
tickpref
box off
xlabel('Q30 of Neuron 1')
ylabel('Q30 of Neuron 2')
title(sprintf('Q30: r = %.2f',r_Q30))
set(gca,'xscale', 'log', 'yscale', 'log');
xlim([min(Q(:,5))*0.95 max(Q(:,5))*1.05])
ylim([min(Q(:,6))*0.95 max(Q(:,6))*1.05])

subplot(7,4,[19 20]);
edges = linspace(0,1.5,21);
n = histc(Q_30diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Q30 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_Q30_diff, med_Q30_diff,std_Q30_diff))
xlabel('Octaves')
tickpref
box off
ylim([0 1.05*max(n)])

subplot(7,4,[25 26]);
hold on
plot([0.05 10],[0.05 10], 'k')
scatter(Q(:,7),Q(:,8),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6], 'SizeData', 18)
xlim([0.05 10])
ylim([0.05 10])
hold off
tickpref
box off
xlabel('Q40 of Neuron 1')
ylabel('Q40 of Neuron 2')
title(sprintf('Q40: r = %.2f',r_Q40))
set(gca,'xscale', 'log', 'yscale', 'log');
xlim([min(Q(:,7))*0.95 max(Q(:,7))*1.05])
ylim([min(Q(:,8))*0.95 max(Q(:,8))*1.05])

subplot(7,4,[27 28])
edges = linspace(0,1.5,21);
n = histc(Q_40diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Q40 Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_Q40_diff, med_Q40_diff,std_Q40_diff))
xlabel('Octaves')
tickpref
box off
ylim([0 1.05*max(n)])

print_mfilename(mfilename);



% Latency 

lat_diff = abs(latency(:,2) - latency(:,1));
med_lat_diff = median(lat_diff);
mean_lat_diff = mean(lat_diff);
std_lat_diff = std(lat_diff);
R_lat = corrcoef(latency(:,1),latency(:,2));
r_lat = R_lat(2);

figure
set(gcf,'Position',[100 100 700 700])
subplot(5,2,[1 6]);
plot(0:40, 0:40, 'k')
hold on
scatter(latency(:,1),latency(:,2),'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold off
tickpref
box off
xlabel('Latency of Neuron 1')
ylabel('Latency of Neuron 2')
title(sprintf('r = %.2f',r_lat))
axis([0 40 0 40])
suptitle('Latency')

subplot(5,2,[9 10]);
edges = linspace(0,10,21);
n = histc(lat_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Latency Difference: Mean = %.2f, Median = %.2f, SD = %.2f, n = 724',mean_lat_diff, med_lat_diff,std_lat_diff))
xlabel('Milliseconds')
tickpref
box off
ylim([0 1.05*max(n)])

print_mfilename(mfilename);




function vs_plot_pairs_fra_fr_scatter_hist(fr)

% Firing Rate
fr_diff = abs(fr(:,2) - fr(:,1));
med_fr_diff = median(fr_diff);
mean_fr_diff = mean(fr_diff);
std_fr_diff = std(fr_diff);
R_fr = corrcoef(fr(:,1),fr(:,2));
r_fr = R_fr(2);

figure
% subplot(5,2,[1 6]);
subplot(2,1,1);
hold on
plot([0 100],[0 100],'k','LineWidth',1.5)
[larger, smaller] = vs_largersmaller(fr(:,1), fr(:,2));
R_fr = corrcoef(larger,smaller);
R_fr(2)

% scatter(fr(:,1),fr(:,2),25, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
scatter(larger,smaller,25, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
hold off
tickpref
box off
xlabel('Firing Rate Neuron 1')
ylabel('Firing Rate Neuron 2')
axis([0 100 0 100]);
xtick = 0:20:100;
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
suptitle('Firing Rates')
title(sprintf('r = %.2f',r_fr))

% subplot(5,2,[9 10]);
subplot(2,1,2);
edges = linspace(0,60,20);
n = histc(fr_diff, edges);
n = n./sum(n);
hb = bar(edges,n,'histc');
set(hb, 'FaceColor', [0.6 0.6 0.6])
title(sprintf('Mean = %.2f, Median = %.2f, SD = %.2f, n = %.0f',...
   mean_fr_diff, med_fr_diff,std_fr_diff, length(fr(:,1))))
xlabel('FR Difference (sp/s)')
tickpref
box off
xlim([min(edges) max(edges)]);
ytick = 0:0.1:0.3;
ylim([min(ytick) max(ytick)])
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gcf,'Position',[177   485   321   404]);
print_mfilename(mfilename);

return




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
rand_median

index = find(rand_median < data_median);
pval = length(index) / nreps;


subplot(1,2,2);
hold on;
[n,x] = hist(rand_median,50);
edges = center2edge(x);
edges
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













