function iccp_plot_fiofit_crosscorr_pos_neg(ccpairs)





close all;


% % Get correlation analysis results
% pd_pos = [ccpairs.pd_pos];
% hw_pos = [ccpairs.hw_pos];
% sigfeature_pos = [ccpairs.sigfeature_pos];
% sigfeature_pos = logical(sigfeature_pos);
% ccc_pos = [ccpairs.ccc_pos];
% ccc_pos(ccc_pos < 0) = 0.0001;
% 
% pd_neg = [ccpairs.pd_neg];
% hw_neg = [ccpairs.hw_neg];
% sigfeature_neg = [ccpairs.sigfeature_neg];
% sigfeature_neg = logical(sigfeature_neg);
% ccc_neg = [ccpairs.ccc_neg];
% ccc_neg(ccc_neg < 0) = 0.0001;
% 
% 
% % Possibilities: Positive only peaks, Negative only peaks, Both Positive and Negative peaks
% % Any significant peak = sum of the previous options
% 
% index_pos_only = sigfeature_pos & ~sigfeature_neg;
% index_neg_only = ~sigfeature_pos & sigfeature_neg;
% index_pos_neg = sigfeature_pos & sigfeature_neg;
% index_any = sigfeature_pos | sigfeature_neg;



[pdPos, hwPos, cccPos, pdNeg, hwNeg, cccNeg, sigPos, sigNeg] = ...
    ccpairs_to_sigfeature(ccpairs);

% Possibilities: Positive only peaks, Negative only peaks, Both Positive and Negative peaks
% Any significant peak = sum of the previous options

index_pos_only = sigPos & ~sigNeg;
index_neg_only = ~sigPos & sigNeg;
index_pos_neg = sigPos & sigNeg;
index_any = sigPos | sigNeg;


fprintf('\n');
fprintf('Positive only peaks: %.0f\n', sum(index_pos_only) );
fprintf('Negative only peaks: %.0f\n', sum(index_neg_only) );
fprintf('Positive and Negative peaks: %.0f\n', sum(index_pos_neg) );
fprintf('Any peaks: %.0f\n', sum(index_any) );
fprintf('\n');


% Get nonlinearity metrics
theta1 = [ccpairs.theta1];
theta2 = [ccpairs.theta2];
theta = [theta1(:) theta2(:)];

sigma1 = [ccpairs.sigma1];
sigma2 = [ccpairs.sigma2];
sigma = [sigma1(:) sigma2(:)];

nmse1 = [ccpairs.nmse1];
index_sig = find(nmse1 <= 0.2);
index_nonsig = find(nmse1 > 0.2);
nmse1(index_sig) = 1;
nmse1(index_nonsig) = 0;
nmse1 = logical(nmse1);

nmse2 = [ccpairs.nmse2];
index_sig = find(nmse2 <= 0.2);
index_nonsig = find(nmse2 > 0.2);
nmse2(index_sig) = 1;
nmse2(index_nonsig) = 0;
nmse2 = logical(nmse2);

nmse = nmse1 & nmse2;


r21 = [ccpairs.r21];
index_sig = find(r21 >=0.8);
index_nonsig = find(r21 < 0.8);
r21(index_sig) = 1;
r21(index_nonsig) = 0;
r21 = logical(r21);

r22 = [ccpairs.r22];
index_sig = find(r22 >=0.8);
index_nonsig = find(r22 < 0.8);
r22(index_sig) = 1;
r22(index_nonsig) = 0;
r22 = logical(r22);

r2 = r21 & r22;

nmse_pos = nmse(index_pos_only);
nmse_neg = nmse(index_neg_only);
nmse_pos_neg = nmse(index_pos_neg);

r2_pos = r2(index_pos_only);
r2_neg = r2(index_neg_only);
r2_pos_neg = r2(index_pos_neg);



% Split Nonlinearity Threshold (theta) into the three correlation groups
theta_pos = theta(index_pos_only,:);
theta_neg = theta(index_neg_only,:);
theta_pos_neg = theta(index_pos_neg,:);

[theta_larger_pos, theta_smaller_pos] = iccp_largersmaller(theta_pos(:,1),theta_pos(:,2));
theta_larger_pos = theta_larger_pos(nmse_pos);
theta_smaller_pos = theta_smaller_pos(nmse_pos);
thetadiff_pos = abs( theta_larger_pos - theta_smaller_pos );

[theta_larger_neg, theta_smaller_neg] = iccp_largersmaller(theta_neg(:,1),theta_neg(:,2));
theta_larger_neg = theta_larger_neg(nmse_neg);
theta_smaller_neg = theta_smaller_neg(nmse_neg);
thetadiff_neg = abs( theta_larger_neg - theta_smaller_neg );

[theta_larger_pos_neg, theta_smaller_pos_neg] = iccp_largersmaller(theta_pos_neg(:,1),theta_pos_neg(:,2));
theta_larger_pos_neg = theta_larger_pos_neg(nmse_pos_neg);
theta_smaller_pos_neg = theta_smaller_pos_neg(nmse_pos_neg);
thetadiff_pos_neg = abs( theta_larger_pos_neg - theta_smaller_pos_neg );







% Make a random distribution, sample, compute paired differences, and plot.
dataUnq = unique(theta);
indexRand = ceil(length(dataUnq) * rand(1,size(theta,1)));
sampleRand1 = dataUnq(indexRand);

indexRand = ceil(length(dataUnq) * rand(1,size(theta,1)));
sampleRand2 = dataUnq(indexRand);

[larger,smaller] = iccp_largersmaller(sampleRand1,sampleRand2);
thetadiffRand = abs( larger - smaller );


pval = ranksum(thetadiffRand, thetadiff_pos)
pval = ranksum(thetadiffRand, thetadiff_pos_neg)


fprintf('\n');
label = 'FIO Theta';
iccp_ccc_pos_neg_summary(thetadiff_pos, thetadiff_neg, thetadiff_pos_neg, label);

fprintf('\n');
fprintf('Sig fit, Positive only peaks: %.0f\n', sum(nmse_pos) );
fprintf('Sig fit, Negative only peaks: %.0f\n', sum(nmse_neg) );
fprintf('Sig fit, Positive and Negative peaks: %.0f\n', sum(nmse_pos_neg) );
fprintf('\n');






% Split Nonlinearity Threshold (theta) into the three correlation groups
sigma_pos = sigma(index_pos_only,:);
sigma_neg = sigma(index_neg_only,:);
sigma_pos_neg = sigma(index_pos_neg,:);

[sigma_larger_pos, sigma_smaller_pos] = iccp_largersmaller(sigma_pos(:,1),sigma_pos(:,2));
sigma_larger_pos = sigma_larger_pos(nmse_pos);
sigma_smaller_pos = sigma_smaller_pos(nmse_pos);
sigmadiff_pos = abs( sigma_larger_pos - sigma_smaller_pos );

[sigma_larger_neg, sigma_smaller_neg] = iccp_largersmaller(sigma_neg(:,1),sigma_neg(:,2));
sigma_larger_neg = sigma_larger_neg(nmse_neg);
sigma_smaller_neg = sigma_smaller_neg(nmse_neg);
sigmadiff_neg = abs( sigma_larger_neg - sigma_smaller_neg );

[sigma_larger_pos_neg, sigma_smaller_pos_neg] = iccp_largersmaller(sigma_pos_neg(:,1),sigma_pos_neg(:,2));
sigma_larger_pos_neg = sigma_larger_pos_neg(nmse_pos_neg);
sigma_smaller_pos_neg = sigma_smaller_pos_neg(nmse_pos_neg);
sigmadiff_pos_neg = abs( sigma_larger_pos_neg - sigma_smaller_pos_neg );


% Randomly order the values
[absc_theta_pos, ord_theta_pos] = ...
    iccp_randomize_columns(theta_pos(:,1), theta_pos(:,2));

[absc_theta_pos_neg, ord_theta_pos_neg] = ...
    iccp_randomize_columns(theta_pos_neg(:,1), theta_pos_neg(:,2));



[absc_sigma_pos, ord_sigma_pos] = ...
    iccp_randomize_columns(sigma_pos(:,1), sigma_pos(:,2));

[absc_sigma_pos_neg, ord_sigma_pos_neg] = ...
    iccp_randomize_columns(sigma_pos_neg(:,1), sigma_pos_neg(:,2));











% Make a random distribution, sample, compute paired differences, and plot.
dataUnq = unique(sigma);
indexRand = ceil(length(dataUnq) * rand(1,size(theta,1)));
sampleRand1 = dataUnq(indexRand);

indexRand = ceil(length(dataUnq) * rand(1,size(theta,1)));
sampleRand2 = dataUnq(indexRand);

[larger,smaller] = iccp_largersmaller(sampleRand1,sampleRand2);
sigmadiffRand = abs( larger - smaller );

pval = ranksum(sigmadiffRand, sigmadiff_pos)
pval = ranksum(sigmadiffRand, sigmadiff_pos_neg)



fprintf('\n');
label = 'FIO Sigma';
iccp_ccc_pos_neg_summary(sigmadiff_pos, sigmadiff_neg, sigmadiff_pos_neg, label);




cmap = brewmaps('blues', 4);
cmap = cmap(1:3,:);
markersize = 4;


figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theta Pairwise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1);
hold on;
xytick = -1:4;
xlim([min(xytick) max(xytick)]);
ylim([min(xytick) max(xytick)]);
plot(xlim, ylim, 'k-');


plot(absc_theta_pos, ord_theta_pos, 'o', ...
    'color', 'k', ...
    'markersize', markersize, ...
    'markeredgecolor', 'k');

gray_face = 0.6 * ones(1,3);
gray_edge = 0.3 * ones(1,3);

plot(absc_theta_pos_neg, ord_theta_pos_neg, 's', ...
    'color', gray_face, ...
    'markersize', markersize, ...
    'markerfacecolor', gray_face, ...
    'markeredgecolor', gray_edge);




tickpref;
xlabel('Theta (SD) Neuron 1');
ylabel('Theta (SD) Neuron 2');
set(gca,'xtick', xytick, 'xticklabel', xytick);
set(gca,'ytick', xytick, 'yticklabel', xytick);
subplot_label(gca,'A');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theta Difference Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2);
hold on;
edges = 0:0.5:4;

n = histc(thetadiff_pos, edges);
pdf_pos = n ./ sum(n);
hp = plot(edges, pdf_pos, 'ko-', ...
    'markersize', markersize, ...
    'markerfacecolor', 'k');
set(hp, 'color', 'k');

n = histc(thetadiff_pos_neg, edges);
pdf_pos_neg = n ./ sum(n);
hp = plot(edges, pdf_pos_neg, 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', 0.6*ones(1,3), ...
    'markeredgecolor', 0.3*ones(1,3));
set(hp, 'color', 0.6*ones(1,3));




n = histc(thetadiffRand, edges);
pdfRand = n ./ sum(n);
hp = plot(edges, pdfRand, 'k-');
set(hp, 'color', 'k', 'linewidth', 2);








hold on;
box off;
tickpref;
ydiffmax = [0.5];
ydifftick = linspace(0,1,11); 
xdifftick = edges(1:2:end);
set(gca,'xtick', xdifftick, 'xticklabel', xdifftick);
set(gca,'ytick', ydifftick, 'yticklabel', ydifftick);
range = max(edges) - min(edges);
xlim([min(edges)-0.025*range max(edges)+0.025*range]);
ylim([0 0.4]);
ylabel('Proportion');
xlabel('Theta Difference (SD)');
legend('EP only', 'EP+ST','Rand');
subplot_label(gca,'B');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma Pairwise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3);
hold on;
xytick = 0.01:5;
xytick = [0.01 0.1 1 5];
xlim([min(xytick) max(xytick)]);
ylim([min(xytick) max(xytick)]);
plot(xlim, ylim, 'k-');


p = 1;
plot(absc_sigma_pos.^(p), ord_sigma_pos.^(p), 'o', ...
    'color', 'k', ...
    'markersize', markersize, ...
    'markeredgecolor', 'k');

gray_face = 0.6 * ones(1,3);
gray_edge = 0.3 * ones(1,3);

plot(absc_sigma_pos_neg.^(p), ord_sigma_pos_neg.^(p), 's', ...
    'color', gray_face, ...
    'markersize', markersize, ...
    'markerfacecolor', gray_face, ...
    'markeredgecolor', gray_edge);

set(gca,'xscale', 'log', 'yscale', 'log');




tickpref;
xlabel('Sigma (SD) Neuron 1');
ylabel('Sigma (SD) Neuron 2');
set(gca,'xtick', xytick, 'xticklabel', xytick);
set(gca,'ytick', xytick, 'yticklabel', xytick);
subplot_label(gca,'C');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma Difference Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(2,2,4);
hold on;
edges = 0:0.5:4;


n = histc(sigmadiff_pos, edges);
pdf_pos = n ./ sum(n);
hp = plot(edges, pdf_pos, 'ko-', ...
    'markersize', markersize, ...
    'markerfacecolor', 'k');
set(hp, 'color', 'k');

n = histc(sigmadiff_pos_neg, edges);
pdf_pos_neg = n ./ sum(n);
hp = plot(edges, pdf_pos_neg, 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', 0.6*ones(1,3), ...
    'markeredgecolor', 0.3*ones(1,3));
set(hp, 'color', 0.6*ones(1,3));


n = histc(sigmadiffRand, edges);
pdfRand = n ./ sum(n);
hp = plot(edges, pdfRand, 'k-');
set(hp, 'color', 'k', 'linewidth', 2);


hold on;
box off;
tickpref;
ydiffmax = [0.6];
ydifftick = linspace(0,1,11); 
xdifftick = edges(1:2:end);
set(gca,'xtick', xdifftick, 'xticklabel', xdifftick);
set(gca,'ytick', ydifftick, 'yticklabel', ydifftick);
range = max(edges) - min(edges);
xlim([min(edges)-0.025*range max(edges)+0.025*range]);
ylim([0 0.5]);
ylabel('Proportion');
xlabel('Sigma Difference (SD)');
subplot_label(gca,'D');


set(gcf,'position', [439 544 438 354]);
print_mfilename(mfilename);


% Estimate Theta correlation values based on random sampling
logtransform = 0;
[rpop_med, rpop_ci, pval] = ...
    iccp_pairwise_corr_rand_test(absc_theta_pos, ord_theta_pos, logtransform);
fprintf('\n');
fprintf(sprintf('%s Pos pairs - randomized\n','Theta'));
fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
fprintf('\n\n');


logtransform = 0;
[rpop_med, rpop_ci, pval] = ...
    iccp_pairwise_corr_rand_test(absc_theta_pos_neg, ord_theta_pos_neg, logtransform);
fprintf('\n');
fprintf(sprintf('%s Pos-Neg pairs - randomized\n','Theta'));
fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
fprintf('\n\n');


% Estimate Sigma correlation values based on random sampling
logtransform = 0;
[rpop_med, rpop_ci, pval] = ...
    iccp_pairwise_corr_rand_test(absc_sigma_pos, ord_sigma_pos, logtransform);
fprintf('\n');
fprintf(sprintf('%s Pos pairs - randomized\n','Sigma'));
fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
fprintf('\n\n');


logtransform = 0;
[rpop_med, rpop_ci, pval] = ...
    iccp_pairwise_corr_rand_test(absc_sigma_pos_neg, ord_sigma_pos_neg, logtransform);
fprintf('\n');
fprintf(sprintf('%s Pos-Neg pairs - randomized\n','Sigma'));
fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
fprintf('\n\n');



return;




sigma_diff = abs( sigma(:,1) - sigma(:,2) );
theta_diff = abs( theta(:,1) - theta(:,2) );
asi_diff = abs( asi(:,1) - asi(:,2) );

edges1 = linspace(floor(min(sigmah)),ceil(max(sigmah)),21);
edges2 = linspace(floor(min(thetah)),ceil(max(thetah)),21);

edges3 = linspace(floor(min(sigma_diff)),ceil(max(sigma_diff)),21);
edges3 = linspace(0,3,21);
edges3 = 0:0.2:4;

edges4 = linspace(floor(min(theta_diff)),ceil(max(theta_diff)),21);
edges4 = linspace(0,4,21);
edges4 = 0:0.2:4;

edges5 = linspace(floor(min(asi_diff)),ceil(max(asi_diff)),21);
edges5 = linspace(0,1,21);
edges5 = 0:0.025:0.5;

label = {'Sigma', 'Theta', 'Asymmetry Index', 'Sigma Difference (SD)', 'Theta Difference (SD)', 'ASI Difference'};

xtick = {edges1, edges2, edges3(1:3:end), edges4(1:3:end), edges5(1:3:end)};

close all;

markersize = 20;


% Histograms of Differences in Sigma, Theta, Asymmetry Index
figure;
subplot(2,2,1)
hold on;
[larger, smaller] = vs_largersmaller(theta(:,1), theta(:,2));
index = larger<4 &larger > -1 & smaller<4 & smaller > -1;
scatter(larger(index), smaller(index), markersize, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel(label{2})
box off
tickpref
xlabel('Theta (SD) Neuron 1');
ylabel('Theta (SD) Neuron 2');
title('Theta')
xytick = -1:4;
set(gca,'xtick', xytick, 'xticklabel', xytick);
set(gca,'ytick', xytick, 'yticklabel', xytick);
xlim([min(xytick) max(xytick)]);
ylim([min(xytick) max(xytick)]);
plot(xlim, ylim, 'k-');

subplot(2,2,2)
N = histc(theta_diff,edges4);
N = N./sum(N);
hb = bar(edges4,N,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
xlabel(label{5})
set(gca,'xtick', 0:4, 'xticklabel', 0:4);
xlim([min(edges4) max(edges4)]);
box off
tickpref
set(gca,'ytick', 0:0.05:0.15, 'yticklabel', 0:0.05:0.15);
ylim([0 0.16]);



fprintf('\n');
fprintf('Theta difference\n');
[r,p] = corrcoef(larger, smaller);
mean_diff = mean(theta_diff);
median_diff = median(theta_diff);
sd_diff = std(theta_diff);
se_diff = std(theta_diff) ./ sqrt(length(theta_diff));
fprintf('Mn=%.4f, Md=%.4f, SD=%.4f, SE=%.4f\n', ...
mean_diff, median_diff, sd_diff, se_diff);
fprintf('r = %.4f, p = %.4f\n', r(2), p(2) );
[pval, data_median] = pairs_permutation_test(larger, smaller);
fprintf('Theta Permutation pval: %.4f\n', pval);
fprintf('\n\n');


subplot(2,2,3)
[larger, smaller] = vs_largersmaller(sigma(:,1), sigma(:,2));
hold on;
scatter(larger(larger<4.5&smaller<4.5), smaller(larger<4.5&smaller<4.5), markersize,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel('Sigma (SD) Neuron 1');
ylabel('Sigma (SD) Neuron 2');
xytick = 0:0.5:4.5;
set(gca,'xtick', xytick, 'xticklabel', xytick);
set(gca,'ytick', xytick, 'yticklabel', xytick);
xlim([min(xytick) max(xytick)]);
ylim([min(xytick) max(xytick)]);
plot(xlim, ylim, 'k-');
box off
tickpref


subplot(2,2,4)
N = histc(sigma_diff,edges3);
N = N./sum(N);
hb = bar(edges3,N,'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
xlabel(label{4})
set(gca,'xtick', 0:4, 'xticklabel', 0:4);
xlim([min(edges3) max(edges3)]);
box off
tickpref
set(gca,'ytick', 0:0.05:0.25, 'yticklabel', 0:0.05:0.25);
ylim([0 0.25]);

fprintf('Sigma difference\n');
[r,p] = corrcoef(larger, smaller);
mean_diff = mean(sigma_diff);
median_diff = median(sigma_diff);
sd_diff = std(sigma_diff);
se_diff = std(sigma_diff) ./ sqrt(length(sigma_diff));
fprintf('Mn=%.4f, Md=%.4f, SD=%.4f, SE=%.4f\n', ...
mean_diff, median_diff, sd_diff, se_diff);
fprintf('r = %.4f, p = %.4f\n', r(2), p(2) );
[pval, data_median] = pairs_permutation_test(larger, smaller);
fprintf('Sigma Permutation pval: %.4f\n', pval);
fprintf('\n\n');






% subplot(3,2,5)
% hold on;
% [larger, smaller] = vs_largersmaller(asi(:,1), asi(:,2));
% scatter(larger, smaller, markersize, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.5 0.5 0.5])
% xlabel(label{3})
% box off
% tickpref
% xlabel('Neuron 1');
% ylabel('Neuron 2');
% title('Asymmetry Index')
% xytick = 0.2:0.2:1;
% set(gca,'xtick', xytick, 'xticklabel', xytick);
% set(gca,'ytick', xytick, 'yticklabel', xytick);
% xlim([min(xytick) max(xytick)]);
% ylim([min(xytick) max(xytick)]);
% plot(xlim, ylim, 'k-');
% 
% subplot(3,2,6)
% mean_asi_diff = mean(asi_diff);
% median_asi_diff = median(asi_diff);
% sd_asi_diff = std(asi_diff);
% se_asi_diff = std(asi_diff) ./ sqrt(length(asi_diff));
% N = histc(asi_diff,edges5);
% N = N./sum(N);
% hb = bar(edges5,N,'histc');
% set(hb, 'facecolor', [0.7 0.7 0.7]);
% xlabel(label{6})
% set(gca,'xtick', 0:0.1:0.5, 'xticklabel', 0:0.1:0.5);
% xlim([0 0.5]);
% title(sprintf('Mn=%.2f,Md=%.2f,SD=%.2f,SE=%.2f', ...
% mean_asi_diff, median_asi_diff, sd_asi_diff, se_asi_diff))
% box off
% tickpref
% ylim([0 0.32]);
print_mfilename(mfilename);
set(gcf,'position', [1125 301 438 410]);


return;



% Sigma Versus Theta Population Plot
figure;
index = thetah > -5 & thetah < 5 & sigmah >= 0 & sigmah < 4.5;
subplot(2,2,3);
scatter(sigmah(index), thetah(index),markersize,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.5 0.5 0.5])
xlabel(label{1})
ylabel(label{2})
box off
tickpref
xlim([0 4.5]);
ylim([-5 5]);
set(gca,'xtick', 0:0.5:4.5, 'xticklabel', 0:0.5:4.5);
set(gca,'ytick', -5:5, 'yticklabel', -5:5);


subplot(2,2,1);
edges = 0:0.5:4.5;
N = histc(sigmah(index), edges);
N = N./sum(N);
hb = bar(edges, N, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
box off
tickpref
xlim([0 max(edges)]);
ylim([0 0.35]);
set(gca,'xtick', edges, 'xticklabel', edges);
set(gca,'ytick', 0:0.05:0.35, 'yticklabel', 0:0.05:0.35);
ylabel('Proportion');


fprintf('\nSigma Population\n');
data = sigmah(index);
fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
   length(data),mean(data), std(data), std(data)./sqrt(length(data)));
fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
   median(data), mad(data), max(data));



subplot(2,2,4);
edges = linspace(-5,5,11);
N = histc(thetah(index), edges);
N = N./sum(N);
hb = barh(edges, N, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
box off
tickpref
xlim([0 0.4]);
ylim([-5 5]);
set(gca,'ytick', -5:5, 'yticklabel', -5:5);
set(gca,'xtick', 0:0.1:0.4, 'xticklabel', 0:0.1:0.4);
xlabel('Proportion');

fprintf('\nTheta Population\n');
data = thetah(index);
fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
   length(data),mean(data), std(data), std(data)./sqrt(length(data)));
fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
   median(data), mad(data), max(data));

set(gcf,'position',[1127         512         438         426]);
print_mfilename(mfilename);



% Plotting Histograms of all Sigma/Theta Values

return;

close all;

figure;
subplot(2,2,[1 2])
meansig = mean(sigmah);
mediansig = median(sigmah);
N = histc(sigmah, edges1);
N = N./sum(N);
hb = bar(edges1, N, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
xlabel(label{1})
set(gca,'xtick', xtick{1}, 'xticklabel', xtick{1});
xlim([min(edges1) max(edges1)]);
title(sprintf('Sigma Population Data, Mean = %.2f, Median = %.2f',meansig,mediansig)) 
box off
tickpref

subplot(2,2,[3 4])
meantheta = mean(thetah);
mediantheta = median(thetah);
M = histc(thetah, edges2);
M = M./sum(M);
hb = bar(edges2, M, 'histc');
set(hb, 'facecolor', [0.7 0.7 0.7]);
xlabel(label{2})
set(gca,'xtick', xtick{2}, 'xticklabel', xtick{2});
xlim([min(edges2) max(edges2)]);
title(sprintf('Theta Population Data, Mean = %.2f, Median = %.2f',meantheta,mediantheta))
box off
tickpref

return;















