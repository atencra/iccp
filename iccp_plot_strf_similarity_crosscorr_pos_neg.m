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

iccp_plot_ccc_strf_similarity_exc_peak(data);

iccp_plot_ccc_strf_similarity_exc_excsup_peak(data);

iccp_plot_ccc_strf_similarity_exc_excsup_sup_peak(data);



% iccp_ccc_strf_similarity_stats(data);


if nargout == 1
    varargout{1} = data;
elseif nargout > 1
    error('There is only 1 output argument.');
end

   
return;







function data = iccp_get_ccc_strf_similarity(ccpairs)

if ( isfield(ccpairs,'pd_pos') && isfield(ccpairs,'pd_neg') )

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

end


if ( ~isfield(ccpairs,'pd_pos') && isfield(ccpairs,'peakdelay') )

    pd_pos = [ccpairs.peakdelay];
    hw_pos = [ccpairs.halfwidth];
    sigfeature_pos = [ccpairs.significant];
    sigfeature_pos = logical(sigfeature_pos);
    ccc_pos = [ccpairs.ccc];
    ccc_pos(ccc_pos < 0) = 0.0001;

    pd_neg = zeros(size(sigfeature_pos)); 
    hw_neg = zeros(size(sigfeature_pos));
    ccc_neg = zeros(size(sigfeature_pos));

    sigfeature_neg = false(size(sigfeature_pos));

    si = [ccpairs.strfsimilarity];

end

data.sigPos = sigfeature_pos(:);
data.sigNeg = sigfeature_neg(:);

data.pdPos = pd_pos(:);
data.hwPos = hw_pos(:);
data.cccPos = ccc_pos(:);

data.pdNeg = pd_neg(:);
data.hwNeg = hw_neg(:);
data.cccNeg = ccc_neg(:);

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



