function iccpairs_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
% iccpairs_plot_strf_similarity_crosscorr_pos_neg STRF similarity and correlation strength
% 
%    iccpairs_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
% 
%    ccpairs : struct array. Each element is one ICC pair of neurons. Data for
%    STRFs and the correlation parameters for the pair are stored in the
%    fields of each element of the struct array.
% 
%    The function plots a histogram of STRF similarity indices, and then 
%    compares the similarity values to the cross-correlation coefficient 
%    values.
% 
%    Data is broken into multiple groups. The CCC data is broken into two
%    groups: the strength of the positive cross-covariance function peak, and
%    the negative peak. Both strengths are plotted against STRF similarity.
%    Only values for significant peaks are plotted.

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

si = [ccpairs.strfsimilarity];


% Split Nonlinearity similarity into the three correlation groups
si_pos = si(index_pos_only);
si_neg = si(index_neg_only);
si_pos_neg = si(index_pos_neg);


fprintf('\n');
label = 'STRF SI';
iccpairs_ccc_pos_neg_summary(si_pos, si_neg, si_pos_neg, label);


cmap = cschemes('greens', 4);
cmap = cmap(1:3,:);
markersize = 4;
linewidth = 2;
ymax = 0.15;

figure;
subplot(2,1,2);
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
ylabel('Cross Correlation Coefficient');
legend('Facilitative', 'Suppressive');
tickpref;
axis([0 1 0.001 1]);
subplot_label(gca,'B');


cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);
subplot(2,1,1);
hold on;
edges_si = 0:0.05:1;
n = histc(si_pos, edges_si);
pdf_pos = n ./ sum(n);
hp = plot(edges_si(1:end-1), pdf_pos(1:end-1), 'o-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
set(hp, 'color', cmap(1,:));

n = histc(si_neg, edges_si);
pdf_neg = n ./ sum(n);
hp = plot(edges_si(1:end-1), smooth3p(pdf_neg(1:end-1)), 'o-', 'markersize', markersize, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
set(hp, 'color', cmap(2,:));

n = histc(si_pos_neg, edges_si);
pdf_pos_neg = n ./ sum(n);
hp = plot(edges_si(1:end-1), smooth3p(pdf_pos_neg(1:end-1)), 'o-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
set(hp, 'color', cmap(3,:));

box off;
tickpref;
xtick = edges_si(1:2:end);
ytick = 0:0.05:0.15; 
ymax = 1;
xtick = 0:0.2:1;
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
range = max(edges_si) - min(edges_si);
xlim([min(edges_si) max(edges_si)]);
ylim([0 0.16]);
ylabel('Proportion');
xlabel('STRF Similarity');
legend('Pos', 'Neg', 'Pos+Neg');
box off;
tickpref;
subplot_label(gca,'A');

set(gcf,'position', [1181 457 321 454]);
print_mfilename(mfilename);











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
fprintf('\nCCPairs: log10 Pos CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigfeature_pos) )), log10(abs(ccc_pos(sigfeature_pos) )));
fprintf('CCPairs: log10 Pos CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigfeature_pos), log10(ccc_pos(sigfeature_pos) ));
fprintf('CCPairs: log10 Pos CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

fprintf('\n\n');

[r,p] = corrcoef(si(sigfeature_neg), log10(ccc_neg(sigfeature_neg)) );
fprintf('\nCCPairs: log10 Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigfeature_neg) )), log10(abs(ccc_neg(sigfeature_neg) )));
fprintf('CCPairs: log10 Neg CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigfeature_neg), log10(ccc_neg(sigfeature_neg) ));
fprintf('CCPairs: log10 Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));


% xfit = linspace(0.1, 1, 1000);
% [beta,s] = polyfit(ccpairs_si_sig, log10(ccpairs_rho_sig), 1); % linear regression
% yfit = polyval(beta,xfit);
% hp = plot(xfit, 10.^yfit, '-', 'color', cmap(1,:), 'linewidth', linewidth);


   
return;


