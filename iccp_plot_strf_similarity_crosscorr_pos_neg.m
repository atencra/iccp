function varargout = iccp_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
% iccp_plot_strf_similarity_crosscorr_pos_neg STRF similarity and correlation strength
% 
%    iccp_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
% 
%    ccpairs : struct array. Each element is one ICC pair of neurons. Data for
%    STRFs and the correlation parameters for the pair are stored in the
%    fields of each element of the struct array.
% 
%    The function plots a histogram of STRF similarity indices, and then 
%    compares the similarity values to the cross-correlation coefficient (CCC) 
%    values.CCC values were estimated for positive and negative
%    cross-covariance deflections.
% 
%    Data is broken into multiple groups. The CCC data is broken into two
%    groups: the strength of the positive cross-covariance function peak, and
%    the negative peak. Both strengths are plotted against STRF similarity.
%    Only values for significant peaks are plotted.
%
%    data = iccp_plot_strf_similarity_crosscorr_pos_neg(ccpairs) returns
%    the data that were used to make the plots. data is a struct with the
%    following fields: sigPos, sigNeg, pdPos, hwPos, cccPos,
%    pdPos, hwPos, cccPos, and si. The varibles are Nx1 or Nx2. 
%
%    pd*, hw*, ccc*, and si are the peak delay, halfwidth, cross-correlation 
%    coefficient, and strf similarity for each pair of neurons. sigPos and 
%    sigNeg are vectors of 0's and 1's specifying if the covariance function 
%    had a significant feature for that pair of neurons.



narginchk(1,1);

data = iccp_get_ccc_strf_similarity(ccpairs);

[sigPosOnly, sigPosNeg, sigNegOnly] = iccp_ccpairs_sigfeature_index(ccpairs);

iccp_ccpairs_crosscorr_strf_si_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);


% Plots with only a positive cross-cov peak
%iccp_plot_ccc_strf_similarity_exc_peak(data);


% Plots with a positive and a positive/negative peak - publication figure.
%iccp_plot_ccc_strf_similarity_exc_excsup_peak(data);


% Plots with a positive peak, a positive/negative peak, and a negative peak
%iccp_plot_ccc_strf_similarity_exc_excsup_sup_peak(data);


%iccp_ccc_strf_similarity_stats(data);


if nargout == 1
    varargout{1} = data;
elseif nargout > 1
    error('There is only 1 output argument.');
end

   
return;







function data = iccp_get_ccc_strf_similarity(ccpairs)


[pdPos, hwPos, cccPos, pdNeg, hwNeg, cccNeg, sigPos, sigNeg] = ...
    ccpairs_to_sigfeature(ccpairs);

si = [ccpairs.strfsimilarity];

data.sigPos = sigPos(:);
data.sigNeg = sigNeg(:);

data.pdPos = pdPos(:);
data.hwPos = hwPos(:);
data.cccPos = cccPos(:);

data.pdNeg = pdNeg(:);
data.hwNeg = hwNeg(:);
data.cccNeg = cccNeg(:);

data.si = si(:);


return;




function iccp_plot_ccc_strf_similarity_exc_peak(data)

cccPos = data.cccPos;
cccNeg = data.cccNeg;
si = data.si;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);

siPosSig = si(indexPos);
siNegSig = si(indexNeg);
siPosNegSig = si(indexPosNeg);


markersize = 4;


figure;

subplot(2,1,1);
cmap = zeros(3,3);
hold on;

edges_si = 0:0.05:1;
n = histc(siPosSig, edges_si);
pdf_pos = n ./ sum(n);
hp = plot(edges_si(1:end-1), pdf_pos(1:end-1), 's-', ...
'markersize', markersize, ...
'markerfacecolor', [0.5 0.5 0.5], ...
'markeredgecolor', 'k' );
set(hp, 'color', cmap(1,:));
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
legend('EP');
box off;
tickpref;
subplot_label(gca,'A');

ht = title(mfilename);
set(ht, 'interpreter', 'none');




subplot(2,1,2);

cmap = zeros(3,3);

hold on;
plot(siPosSig, cccPosSig, 's', ...
'color', cmap(1,:), ...
'markersize', markersize, ...
'markerfacecolor', [0.5 0.5 0.5], ...
'markeredgecolor', 'k');

ytick = [0.0001 0.001 0.01 0.1 1];
xtick = [0:0.2:1];
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'yscale', 'log');
xlabel('STRF Similarity');
ylabel('Cross Correlation Coefficient');
legend('EP');
tickpref;
axis([0 1 0.001 1]);
subplot_label(gca,'B');


set(gcf,'position', [1181 457 321 454]);

return;




function iccp_plot_ccc_strf_similarity_exc_excsup_peak(data)

cccPos = data.cccPos;
cccNeg = data.cccNeg;
si = data.si;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);

siPosSig = si(indexPos);
siNegSig = si(indexNeg);
siPosNegSig = si(indexPosNeg);


markersize = 4;



if ( sum(indexNeg) )

    figure;

    subplot(2,1,1);

    cmap = brewmaps('blues', 4);
    cmap = cmap(1:3,:);

    hold on;

    edges_si = 0:0.05:1;
    n = histc(siPosSig, edges_si);
    pdf_pos = n ./ sum(n);
    hp = plot(edges_si(1:end-1), pdf_pos(1:end-1), 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', cmap(1,:) );
    set(hp, 'color', cmap(1,:));

    n = histc(siPosNegSig, edges_si);
    pdf_pos_neg = n ./ sum(n);
    hp = plot(edges_si(1:end-1), smooth3p(pdf_pos_neg(1:end-1)), 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', cmap(3,:), ...
    'markeredgecolor', cmap(3,:) );
    set(hp, 'color', cmap(3,:));
    edges_si = 0:0.05:1;


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
    legend('EP', 'EP+SP');
    box off;
    tickpref;
    subplot_label(gca,'A');

    ht = title(mfilename);
    set(ht, 'interpreter', 'none');




    subplot(2,1,2);

    cmap = brewmaps('greens', 4);
    cmap = cmap(1:3,:);


    hold on;

    plot(siPosSig, cccPosSig, 's', ...
        'color', cmap(1,:), ...
        'markersize', markersize, ...
        'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

    plot(siNegSig, cccNegSig, 's', ...
        'color', cmap(3,:), ...
        'markersize', markersize, ...
        'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

    ytick = [0.0001 0.001 0.01 0.1 1];
    xtick = [0:0.2:1];
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    set(gca,'ytick', ytick, 'yticklabel', ytick);
    set(gca,'yscale', 'log');
    xlabel('STRF Similarity');
    ylabel('Cross Correlation Coefficient');
    legend('EP', 'SP');
    tickpref;
    axis([0 1 0.001 1]);
    subplot_label(gca,'B');

    set(gcf,'position', [1181 457 321 454]);

end

return;




function iccp_plot_ccc_strf_similarity_exc_excsup_sup_peak(data)

cccPos = data.cccPos;
cccNeg = data.cccNeg;
si = data.si;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);

siPosSig = si(indexPos);
siNegSig = si(indexNeg);
siPosNegSig = si(indexPosNeg);

markersize = 4;


if ( sum(indexNeg) )

    figure;

    subplot(2,1,1);

    cmap = brewmaps('blues', 4);
    cmap = cmap(1:3,:);

    hold on;

    edges_si = 0:0.05:1;
    n = histc(siPosSig, edges_si);
    pdf_pos = n ./ sum(n);
    hp = plot(edges_si(1:end-1), pdf_pos(1:end-1), 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', cmap(1,:) );
    set(hp, 'color', cmap(1,:));

    n = histc(siNegSig, edges_si);
    pdf_neg = n ./ sum(n);
    hp = plot(edges_si(1:end-1), smooth3p(pdf_neg(1:end-1)), 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', cmap(2,:), ...
    'markeredgecolor', cmap(2,:) );
    set(hp, 'color', cmap(2,:));

    n = histc(siPosNegSig, edges_si);
    pdf_pos_neg = n ./ sum(n);
    hp = plot(edges_si(1:end-1), smooth3p(pdf_pos_neg(1:end-1)), 's-', ...
    'markersize', markersize, ...
    'markerfacecolor', cmap(3,:), ...
    'markeredgecolor', cmap(3,:) );
    set(hp, 'color', cmap(3,:));
    edges_si = 0:0.05:1;


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
    legend('EP', 'SP','EP+SP');
    box off;
    tickpref;
    subplot_label(gca,'A');

    ht = title(mfilename);
    set(ht, 'interpreter', 'none');


    subplot(2,1,2);

    cmap = brewmaps('greens', 4);
    cmap = cmap(1:3,:);

    hold on;

    plot(siPosSig, cccPosSig, 's', ...
    'color', cmap(1,:), ...
    'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', cmap(1,:));

    plot(siNegSig, cccNegSig, 's', ...
    'color', cmap(3,:), ...
    'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
    'markeredgecolor', cmap(3,:));

    ytick = [0.0001 0.001 0.01 0.1 1];
    xtick = [0:0.2:1];
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    set(gca,'ytick', ytick, 'yticklabel', ytick);
    set(gca,'yscale', 'log');
    xlabel('STRF Similarity');
    ylabel('Cross Correlation Coefficient');
    legend('EP', 'SP');
    tickpref;
    axis([0 1 0.001 1]);
    subplot_label(gca,'B');

    set(gcf,'position', [1181 457 321 454]);

end

return;




function iccp_ccc_strf_similarity_stats(data)

si = data.si;
sigPos = data.sigPos;
sigNeg = data.sigNeg;
cccPos = data.cccPos;
cccNeg = data.cccNeg;

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


[r,p] = corrcoef(si(sigPos), log10(cccPos(sigPos)) );
fprintf('\nCCPairs: log10 Pos CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigPos) )), log10(abs(cccPos(sigPos) )));
fprintf('CCPairs: log10 Pos CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigPos), log10(cccPos(sigPos) ));
fprintf('CCPairs: log10 Pos CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

fprintf('\n\n');

[r,p] = corrcoef(si(sigNeg), log10(cccNeg(sigNeg)) );
fprintf('\nCCPairs: log10 Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigNeg) )), log10(abs(cccNeg(sigNeg) )));
fprintf('CCPairs: log10 Neg CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigNeg), log10(cccNeg(sigNeg) ));
fprintf('CCPairs: log10 Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));


% xfit = linspace(0.1, 1, 1000);
% [beta,s] = polyfit(ccpairs_si_sig, log10(ccpairs_rho_sig), 1); % linear regression
% yfit = polyval(beta,xfit);
% hp = plot(xfit, 10.^yfit, '-', 'color', cmap(1,:), 'linewidth', linewidth);


return;






function [sigPosOnly, sigPosNeg, sigNegOnly] = iccp_ccpairs_sigfeature_index(ccpairs)


pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);

sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);


if ( isfield(ccpairs,'pd_pos') && isfield(ccpairs,'pd_neg') )
    pdNeg = [ccpairs.pd_neg];
    pdNeg = abs(pdNeg);

    sigNeg = [ccpairs.sigfeature_neg];

    sigPosOnly = sigPos & ~sigNeg;

    sigNeg0 = logical(sigNeg & pdNeg > 0); % for positive/negative cross-cov
    sigNeg15 = logical(sigNeg & pdNeg > 1.5); % for negative only cross-cov

    sigNegOnly = ~sigPos & sigNeg15;
    sigPosNeg = sigPos & sigNeg0;

    fprintf('\n');
    fprintf('Positive only peaks: %.0f\n', sum(sigPosOnly) );
    fprintf('Negative only peaks: %.0f\n', sum(sigNegOnly) );
    fprintf('Positive and Negative peaks: %.0f\n', sum(sigPosNeg) );
    fprintf('Any peaks: %.0f\n', sum([sigPosOnly(:); sigNegOnly(:); sigPosNeg(:)]) );
    fprintf('\n');
end


if ( ~isfield(ccpairs,'pd_pos') && isfield(ccpairs,'peakdelay') )
    sigPosOnly = sigPos;
    sigNegOnly = false(size(sigPos));
    sigPosNeg = false(size(sigPosNeg));
end

return;




function iccp_ccpairs_crosscorr_strf_si_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)


% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);


% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);

hwPos = [ccpairs.hw_pos];

cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;


% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);

hwNeg = [ccpairs.hw_neg];

cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;


si = [ccpairs.strfsimilarity];

[length(cccPos) length(cccNeg) length(si)]


fprintf('\n');
fprintf('STRF Similarity population\n');
fprintf('N = %.0f\n', length(si));
fprintf('MN = %.4f\n', mean(si(~isnan(si))));
fprintf('SD = %.4f\n', std(si(~isnan(si) ) ) );
fprintf('MD = %.4f\n', median(si(~isnan(si))));
fprintf('MAD = %.4f\n', madstat(si(~isnan(si))));
fprintf('MX = %.4f\n', max(si(~isnan(si))) );


fprintf('\n');
fprintf('STRF Similarity: Pos Only\n');
siPos = si(sigPosOnly);
fprintf('N = %.0f\n', length(siPos));
fprintf('MN = %.4f\n', mean(siPos));
fprintf('SD = %.4f\n', std(siPos));
fprintf('MD = %.4f\n', median(siPos));
fprintf('MAD = %.4f\n', madstat(siPos));
fprintf('MX = %.4f\n', max(siPos));


fprintf('\n');
fprintf('STRF Similarity: Pos Neg Only\n');
siPosNeg = si(sigPosNeg);
fprintf('N = %.0f\n', length(siPosNeg));
fprintf('MN = %.4f\n', mean(siPosNeg));
fprintf('SD = %.4f\n', std(siPosNeg));
fprintf('MD = %.4f\n', median(siPosNeg));
fprintf('MAD = %.4f\n', madstat(siPosNeg));
fprintf('MX = %.4f\n', max(siPosNeg));


fprintf('\n');
fprintf('STRF Similarity: Neg Only\n');
siNeg = si(sigNegOnly);
fprintf('N = %.0f\n', length(siNeg));
fprintf('MN = %.4f\n', mean(siNeg));
fprintf('SD = %.4f\n', std(siNeg));
fprintf('MD = %.4f\n', median(siNeg));
fprintf('MAD = %.4f\n', madstat(siNeg));
fprintf('MX = %.4f\n', max(siNeg));



pPos_vs_PosNeg = ranksum(siPos, siPosNeg);
pPos_vs_Neg = ranksum(siPos, siNeg);
pPosNeg_vs_Neg = ranksum(siPosNeg, siNeg);

fprintf('\n');
fprintf('SI: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('SI: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('SI: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


fprintf('\n');

[r,p] = corrcoef(si(sigPosOnly), log10(cccPos(sigPosOnly)) );
fprintf('\nCCPairs: log10 Pos Only CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigPosOnly) )), log10(abs(cccPos(sigPosOnly) )));
fprintf('CCPairs: log10 Pos Only CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigPosOnly), log10(cccPos(sigPosOnly) ));
fprintf('CCPairs: log10 Pos Only CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

fprintf('\n');

[r,p] = corrcoef(si(sigPosNeg), log10(cccNeg(sigPosNeg)) );
fprintf('\nCCPairs: log10 Pos Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigPosNeg) )), log10(abs(cccNeg(sigPosNeg) )));
fprintf('CCPairs: log10 Pos Neg CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigPosNeg), log10(cccNeg(sigPosNeg) ));
fprintf('CCPairs: log10 Pos Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));


fprintf('\n\n');







% cccPos = data.cccPos;
% cccNeg = data.cccNeg;
% si = data.si;
% 
% indexPos = data.sigPos;
% indexNeg = data.sigNeg;
% indexPosNeg = indexPos & indexNeg;
% 
% cccPosSig = cccPos(indexPos);
% cccNegSig = cccNeg(indexNeg);
% 
% siPosSig = si(indexPos);
% siNegSig = si(indexNeg);
% siPosNegSig = si(indexPosNeg);
% 
% 
% markersize = 4;
% 
% 
% 
% if ( sum(indexNeg) )
% 
%     figure;
% 
%     subplot(2,1,1);
% 
%     cmap = brewmaps('blues', 4);
%     cmap = cmap(1:3,:);
% 
%     hold on;
% 
%     edges_si = 0:0.05:1;
%     n = histc(siPosSig, edges_si);
%     pdf_pos = n ./ sum(n);
%     hp = plot(edges_si(1:end-1), pdf_pos(1:end-1), 's-', ...
%     'markersize', markersize, ...
%     'markerfacecolor', cmap(1,:), ...
%     'markeredgecolor', cmap(1,:) );
%     set(hp, 'color', cmap(1,:));
% 
%     n = histc(siPosNegSig, edges_si);
%     pdf_pos_neg = n ./ sum(n);
%     hp = plot(edges_si(1:end-1), smooth3p(pdf_pos_neg(1:end-1)), 's-', ...
%     'markersize', markersize, ...
%     'markerfacecolor', cmap(3,:), ...
%     'markeredgecolor', cmap(3,:) );
%     set(hp, 'color', cmap(3,:));
%     edges_si = 0:0.05:1;
% 
% 
%     box off;
%     tickpref;
%     xtick = edges_si(1:2:end);
%     ytick = 0:0.05:0.15; 
%     ymax = 1;
%     xtick = 0:0.2:1;
%     set(gca,'xtick', xtick, 'xticklabel', xtick);
%     set(gca,'ytick', ytick, 'yticklabel', ytick);
%     range = max(edges_si) - min(edges_si);
%     xlim([min(edges_si) max(edges_si)]);
%     ylim([0 0.16]);
%     ylabel('Proportion');
%     xlabel('STRF Similarity');
%     legend('EP', 'EP+SP');
%     box off;
%     tickpref;
%     subplot_label(gca,'A');
% 
%     ht = title(mfilename);
%     set(ht, 'interpreter', 'none');
% 
% 
% 
% 
%     subplot(2,1,2);
% 
%     cmap = brewmaps('greens', 4);
%     cmap = cmap(1:3,:);
% 
% 
%     hold on;
% 
%     plot(siPosSig, cccPosSig, 's', ...
%         'color', cmap(1,:), ...
%         'markersize', markersize, ...
%         'markerfacecolor', cmap(1,:), ...
%         'markeredgecolor', cmap(1,:));
% 
%     plot(siNegSig, cccNegSig, 's', ...
%         'color', cmap(3,:), ...
%         'markersize', markersize, ...
%         'markerfacecolor', cmap(3,:), ...
%         'markeredgecolor', cmap(3,:));
% 
%     ytick = [0.0001 0.001 0.01 0.1 1];
%     xtick = [0:0.2:1];
%     set(gca,'xtick', xtick, 'xticklabel', xtick);
%     set(gca,'ytick', ytick, 'yticklabel', ytick);
%     set(gca,'yscale', 'log');
%     xlabel('STRF Similarity');
%     ylabel('Cross Correlation Coefficient');
%     legend('EP', 'SP');
%     tickpref;
%     axis([0 1 0.001 1]);
%     subplot_label(gca,'B');
% 
%     set(gcf,'position', [1181 457 321 454]);
% 
% end





return;




