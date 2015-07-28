function [sigmah, thetah, sigma, theta, asi] = ...
   vs_nonlinearity_parameter_analysis(sigmah, thetah, sigma, theta, asi)
% nonlinearity_parameter_analysis: Looks at nonlinearity parameters from
%                                  *-sta-fio-fitparams.mat files.

% Loops through all sta-fio-fitparams data, and:
% 1) Performs pairwise analysis on sigma, theta, and asymmetry index values.
% 2) Plots difference histograms for sigma, theta, and asymmetry index
%       pairs, as well as scatterplots of the pairs.
% 3) Plots sigma population values vs. theta population values.
% 4) Plots population histograms for sigma and theta values.


if ( nargin == 0 ) % read in all the data from the files
   d = dir('*-sta-fio-fitparams.mat');
   sigma = [];
   theta = [];
   asi = [];
   sigmah = [];
   thetah = [];

   for n = 1:length(d)
      filename = d(n).name;
      load(filename, 'fioFit');

      chan = [fioFit.chan];
      chan_unique = unique(chan);

      for m = 1:length(fioFit)
         if ( fioFit(m).nmse_sta < 0.1 )
            sigmah = [sigmah; fioFit(m).fitParams_sta(3)];
            thetah = [thetah; fioFit(m).fitParams_sta(2)];
         end
      end

       for i = 1:length(chan_unique)

           index = find( chan_unique(i) == chan );

           if length(index)>1

               cmb = nchoosek(index, 2); % determine all possible pairwise combinations

               [nr, nc] = size(cmb);

               for j = 1:nr

                   index1 = cmb(j,1);
                   index2 = cmb(j,2);

                   if fioFit(index1).nmse_sta < 0.1 && fioFit(index2).nmse_sta < 0.1

                       sigma1 = fioFit(index1).fitParams_sta(3);
                       theta1 = fioFit(index1).fitParams_sta(2);
                       sigma2 = fioFit(index2).fitParams_sta(3);
                       theta2 = fioFit(index2).fitParams_sta(2);
                       asi1 = fioFit(index1).asi_sta;
                       asi2 = fioFit(index2).asi_sta;

                       sigma = [sigma; sigma1 sigma2];
                       theta = [theta; theta1 theta2];
                       asi = [asi; asi1 asi2];

                   end

               end % for j

           end % if

       end % for i

   end % for n

end % (if nargin == 0)


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















