function iccp_plot_nonlinearity_asi_si_crosscorr_pos_neg(ccpairs)
% iccp_plot_nonlinearity_asi_si_crosscorr_pos_neg Compare nonlinearity ASI, SI
% 
%    iccpairs_plot_nonlinearity_asi_si_crosscorr_pos_neg(ccpairs)
% 
%    Compares the nonlinearities for a pair of ICC neurons. Two metrics are used:
%    the asymmetry of the nonlinearity and the similarity, or correlation,
%    between the two nonlinearities.
% 
%    Data are split into three groups based on the correlation function:
%    1) Cross-covariance functions with only Positive peaks
%    2) Cross-covariance functions with only Negative 'peaks'
%    3) Cross-covariance functions with only both Positive and Negative peaks
% 
%    ccpairs : struct array, with each element representing one pair of 
%    neurons. The fields of ccpairs hold the nonlinearitiy parameters
%    and the correlation analysis results.

close all;


[pdPos, hwPos, cccPos, pdNeg, hwNeg, cccNeg, sigPos, sigNeg] = ...
    ccpairs_to_sigfeature(ccpairs);

% Possibilities: Positive only peaks, Negative only peaks, Both Positive and Negative peaks
% Any significant peak = sum of the previous options

index_pos_only = sigPos & ~sigNeg;
index_neg_only = ~sigPos & sigNeg;
index_pos_neg = sigPos & sigNeg;
index_any = sigPos | sigNeg;

fprintf('Positive only peaks: %.0f\n', sum(index_pos_only) );
fprintf('Negative only peaks: %.0f\n', sum(index_neg_only) );
fprintf('Positive and Negative peaks: %.0f\n', sum(index_pos_neg) );
fprintf('Any peaks: %.0f\n', sum(index_any) );


% Get nonlinearity metrics
asi1 = [ccpairs.fioasi1];
asi2 = [ccpairs.fioasi2];
asi = [asi1(:) asi2(:)];
si = [ccpairs.fiosimilarity];


% Split Nonlinearity similarity into the three correlation groups
si_pos = si(index_pos_only);
si_neg = si(index_neg_only);
si_pos_neg = si(index_pos_neg);


% Split Nonlinearity asymmetry into the three correlation groups
asi_pos = asi(index_pos_only,:);
asi_neg = asi(index_neg_only,:);
asi_pos_neg = asi(index_pos_neg,:);


% Order asymmetry based on magnitude, and calculate the difference between
% the values for a pair of neurons
[larger_pos,smaller_pos] = iccp_largersmaller(asi_pos(:,1),asi_pos(:,2));
asidiff_pos = abs( larger_pos - smaller_pos );

[larger_neg,smaller_neg] = iccp_largersmaller(asi_neg(:,1),asi_neg(:,2));
asidiff_neg = abs( larger_neg - smaller_neg );

[larger_pos_neg,smaller_pos_neg] = iccp_largersmaller(asi_pos_neg(:,1),asi_pos_neg(:,2));
asidiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );

fprintf('\n');
label = 'FIO ASI';
iccp_ccc_pos_neg_summary(asidiff_pos, asidiff_neg, asidiff_pos_neg, label);



fprintf('\n');
label = 'FIO SI';
iccp_ccc_pos_neg_summary(si_pos, si_neg, si_pos_neg, label);



% Plotting specs
cmap = brewmaps('blues', 4);
cmap = cmap(1:3,:);
markersize = 4;


figure;

subplot(3,1,1);
hold on;
xlim([0.5 1]);
ylim([xlim]);
plot([xlim], [xlim], 'k-');

plot(larger_pos, smaller_pos, 's', 'color', cmap(1,:), ...
   'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
   'markeredgecolor', cmap(1,:));

% plot(larger_neg, smaller_neg, 'o', 'color', cmap(2,:), ...
%    'markersize', markersize, 'markerfacecolor', cmap(2,:), ...
%    'markeredgecolor', cmap(2,:));

plot(larger_pos_neg, smaller_pos_neg, 's', 'color', cmap(3,:), ...
   'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
   'markeredgecolor', cmap(3,:));

tickpref;
tick = 0.5:0.1:1;
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
xlabel('ASI Neuron 1');
ylabel('ASI Neuron 2');
subplot_label(gca,'A');


subplot(3,1,2);
hold on;
edges_asi = 0:0.025:0.5;
n = histc(asidiff_pos, edges_asi);
pdf_pos = n ./ sum(n);
hp = plot(edges_asi, pdf_pos, 's-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
set(hp, 'color', cmap(1,:));

% n = histc(asidiff_neg, edges_asi);
% pdf_neg = n ./ sum(n);
% hp = plot(edges_asi, pdf_neg, 'o-', 'markersize', markersize, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
% set(hp, 'color', cmap(2,:));

n = histc(asidiff_pos_neg, edges_asi);
pdf_pos_neg = n ./ sum(n);
hp = plot(edges_asi, pdf_pos_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
set(hp, 'color', cmap(3,:));

hold on;
box off;
tickpref;
ydiffmax = [0.5];
ydifftick = linspace(0,1,11); 
xdifftick = edges_asi(1:2:end);
set(gca,'xtick', 0:0.1:0.5, 'xticklabel', 0:0.1:0.5);
set(gca,'ytick', 0:0.1:0.4, 'yticklabel', 0:0.1:0.4);
range = max(edges_asi) - min(edges_asi);
xlim([min(edges_asi)-0.025*range max(edges_asi)+0.025*range]);
ylim([0 0.4]);
ylabel('Proportion');
xlabel('ASI Difference');
legend('Pos', 'Pos+Neg');
subplot_label(gca,'B');




subplot(3,1,3);
hold on;
edges_si = 0.4:0.025:1;
n = histc(si_pos, edges_si);
pdf_pos = n ./ sum(n);
hp = plot(edges_si(1:end-1), pdf_pos(1:end-1), 's-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
set(hp, 'color', cmap(1,:));

% n = histc(si_neg, edges_si);
% pdf_neg = n ./ sum(n);
% hp = plot(edges_si(1:end-1), pdf_neg(1:end-1), 'o-', 'markersize', markersize, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
% set(hp, 'color', cmap(2,:));

n = histc(si_pos_neg, edges_si);
pdf_pos_neg = n ./ sum(n);
hp = plot(edges_si(1:end-1), pdf_pos_neg(1:end-1), 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
set(hp, 'color', cmap(3,:));

% [edges_si(:) pdf_pos(:) pdf_neg(:) pdf_pos_neg(:)]

hold on;
box off;
tickpref;
xtick = edges_si(1:2:end);
ytick = linspace(0,1,11); 
ymax = 1;
set(gca,'xtick', 0.4:0.1:1, 'xticklabel', 0.4:0.1:1);
set(gca,'ytick', 0:0.2:0.6, 'yticklabel', 0:0.2:0.6);
range = max(edges_si) - min(edges_si);
xlim([min(edges_si) max(edges_si)]);
ylim([0 0.6]);
ylabel('Proportion');
xlabel('FIO SI');
subplot_label(gca,'C');

% legend('Pos', 'Neg', 'Pos+Neg');


% % Permutation test
% [pval, data_median] = pairs_permutation_test(larger, smaller);
% fprintf('\nAsymmetry Index Difference\n');
% fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
%    length(datadiff),mean(datadiff), std(datadiff), std(datadiff)./sqrt(length(datadiff)));
% fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
%    median(datadiff), mad(datadiff), max(datadiff));
% fprintf('ASI Permutation pval: %.4f\n', pval);
% fprintf('\n\n');




% fprintf('\nFIO Similarity\n');
% fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
%    length(ccpairs_si),mean(ccpairs_si), std(ccpairs_si), std(ccpairs_si)./sqrt(length(ccpairs_si)));
% fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
%    median(ccpairs_si), mad(ccpairs_si), max(ccpairs_si));


set(gcf,'position', [465 207 321 695]);
print_mfilename(mfilename);



return;










