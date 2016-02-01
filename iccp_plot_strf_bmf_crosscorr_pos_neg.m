function varargout = iccp_plot_strf_bmf_crosscorr_pos_neg(ccpairs)
% iccpairs_plot_strf_bmf_crosscorr_pos_neg BMF for significant paired correlations
% 
%    iccp_plot_strf_bmf_crosscorr_pos_neg(ccpairs)
% 
%    Plots best temporal and spectral modulation frequencies for significant
%    cross-covariance functions. The data are split into three groups:
%    1) When cross-covariance functions have only a positive peak
%    2) For only a negative peak
%    3) When there is both a positive and a negative peak
% 
%    The figure shows a plot of the bmfs for each pair of neurons plotted
%    against each other. Then the difference between the bmfs is shown.
%    Last, the bmf difference is plotted against the cross-correlation
%    coefficient for the positive and negative peaks in the cross-covariance
%    function.
% 
%    ccpairs : struct array. Each element is one ICC pair of neurons. Data for
%    STRFs and the correlation parameters for the pair are stored in the
%    fields of each element of the struct array.
%
%    data = iccp_plot_strf_bmf_crosscorr_pos_neg(ccpairs) returns
%    the data that were used to make the plots. data is a struct with the
%    following fields: 'sigPos', 'sigNeg', 'pdPos', 'hwPos', 'cccPos', 
%    'pdNeg', 'hwNeg', 'cccNeg', 'btmf', 'bsmf'.
%
%    pd*, hw*, and ccc* are the peak delay, halfwidth, cross-correlation 
%    coefficient. btmf and bsmf are the best temporal and spectral modulation
%    for each neuron. They are Nx2 arrays, with the data for one
%    pair on each row. sigPos and sigNeg are vectors of 0's and 1's 
%    specifying if the covariance function had a significant feature for 
%    that pair of neurons.




close all;

data = iccp_ccpairs2bmf(ccpairs);

%iccp_plot_bmf_ccc_exc_peak(data);

iccp_plot_bmf_ccc_exc_excsup_peak(data);

%iccp_plot_bmf_ccc_exc_excsup_sup_peak(data);


if ( nargout == 1 )
    varargout{1} = data;
end


return







function data = iccp_ccpairs2bmf(ccpairs)


[pdPos, hwPos, cccPos, pdNeg, hwNeg, cccNeg, sigPos, sigNeg] = ...
    ccpairs_to_sigfeature(ccpairs);

data.sigPos = sigPos(:);
data.sigNeg = sigNeg(:);

data.pdPos = pdPos(:);
data.hwPos = hwPos(:);
data.cccPos = cccPos(:);

data.pdNeg = pdNeg(:);
data.hwNeg = hwNeg(:);
data.cccNeg = cccNeg(:);



% Paired modulation data
btmf1 = [ccpairs.btmf1];
btmf2 = [ccpairs.btmf2];

bsmf1 = [ccpairs.bsmf1];
bsmf2 = [ccpairs.bsmf2];

btmf = [btmf1(:) btmf2(:)];
bsmf = [bsmf1(:) bsmf2(:)];

data.btmf = btmf;
data.bsmf = bsmf;

return;


function iccp_plot_bmf_ccc_exc_peak(data)


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);


ccpairs_btmf = data.btmf;
ccpairs_bsmf = data.bsmf;

spktype_label{1} = 'CCPairs';
spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';
datatype = {'btmf', 'bsmf'};

xaxScatter =  {[12.5 400],[0.05 2.1]};
yaxScatter = xaxScatter;

xscale = {'log', 'log', 'log'};
yscale = {'log', 'log', 'log'};

edges = {linspace(0, 300, 21), linspace(0, 2, 21)};

xtickScatter = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
xtickHist = {[0:60:300], [0:0.4:2]}; 
xtickDiff = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1 2]}; 


xlabelScatter = {'bTMF (Hz)', 'bSMF (cyc/oct)'};
xlabelHist = {'bTMF Difference (Hz)', 'bSMF Difference (cyc/oct)'};
xlabelCCC = xlabelHist;

yaxDiff = {[0 0.4], [0 0.6]};

ytickHist = {0:0.1:0.4, 0:0.2:0.6, linspace(0,0.6,4)}; 

ytickCCC = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1]}; 


xaxCCC = {[0.01 400],[0.001 2]};
yaxCCC = {[0.001 1], [0.001 1]};

markersize = 4;



close all;

for i = 1:length(datatype)


    figure;
    dataspktype = eval(['ccpairs_' datatype{i}]);
    dataspktype_pos = dataspktype(indexPos,:);
    dataspktype_neg = dataspktype(indexNeg,:);
    dataspktype_pos_neg = dataspktype(indexPosNeg,:);

    [larger,smaller] = iccp_largersmaller(dataspktype(:,1),dataspktype(:,2));
    datadiff = abs( larger - smaller );


    [larger_pos,smaller_pos] = ...
        iccp_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));
    datadiff_pos = abs( larger_pos - smaller_pos );

    if ( sum(indexNeg) )
        [larger_neg,smaller_neg] = ...
            iccp_largersmaller(dataspktype_neg(:,1),dataspktype_neg(:,2));
        datadiff_neg = abs( larger_neg - smaller_neg );

        [larger_pos_neg,smaller_pos_neg] = ...
            iccp_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));
        datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );
        iccp_ccc_pos_neg_summary(datadiff_pos, datadiff_neg, datadiff_pos_neg, datatype{i});
    end



    subplot(3,1,1);
    hold on;
    plot([xaxScatter{i}], [yaxScatter{i}], 'k-');


    if ( ~sum(indexNeg) )
        plot(larger_pos, smaller_pos, 's', ...
        'markersize', markersize, ... %'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', 0.25*ones(1,3)); %cmap(1,:));

    else

        cmap = brewmaps('blues', 4);
        cmap = cmap(1:3,:);

        plot(larger_pos, smaller_pos, 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        plot(larger_neg, smaller_neg, 's', 'color', cmap(2,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(2,:), ...
        'markeredgecolor', cmap(2,:));

        plot(larger_pos_neg, smaller_pos_neg, 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));
    end


    set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
    tickpref;
    set(gca,'xtick', xtickScatter{i}, 'xticklabel', xtickScatter{i});
    set(gca,'ytick', xtickScatter{i}, 'yticklabel', xtickScatter{i});

    xlim([xaxScatter{i}]);
    ylim([yaxScatter{i}]);
    xlabel(sprintf('%s %s', xlabelScatter{i}, spktype_xlabel{1}));
    ylabel(sprintf('%s %s', xlabelScatter{i}, spktype_ylabel{1}));
    subplot_label(gca,'A');




    subplot(3,1,2);
    hold on;

    if ( ~sum(indexNeg) )
        n = histc(datadiff_pos, edges{i});
        pdf_pos = n ./ sum(n);
        plot(edges{i}, pdf_pos, 'ks-', 'markersize', 5, ...%markersize, ...
            'markerfacecolor', 0.5*ones(1,3), ...
            'markeredgecolor', zeros(1,3)); 
    else % if ( sum(indexNeg) )
        n = histc(datadiff_pos, edges{i});
        pdf_pos = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos, 's-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
        set(hp, 'color', cmap(1,:));

        n = histc(datadiff_neg, edges{i});
        pdf_neg = n ./ sum(n);
        hp = plot(edges{i}, pdf_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
        set(hp, 'color', cmap(2,:));

        n = histc(datadiff_pos_neg, edges{i});
        pdf_pos_neg = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
        set(hp, 'color', cmap(3,:));
    end

    hold on;
    box off;
    tickpref;
    set(gca,'xtick', xtickHist{i}, 'xticklabel', xtickHist{i});
    set(gca,'ytick', ytickHist{i}, 'yticklabel', ytickHist{i});
    range = max(edges{i}) - min(edges{i});
    xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
    ylim([yaxDiff{i}]);

    xlabel(xlabelHist{i});
    ylabel('Proportion');

    if ( sum(indexNeg) )
    legend('EP', 'SP','EP+SP');
    else
    legend('EP');
    end
    subplot_label(gca,'B');



    subplot(3,1,3);
    hold on;

    if ( ~sum(indexNeg) )
        plot(datadiff(indexPos), cccPos(indexPos), 's', ...
        'markersize', markersize, ... %'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', 0.25*ones(1,3)); %cmap(1,:));

        x = log10(datadiff(indexPos));
        y = log10(cccPos(indexPos));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

    else % if ( sum(indexNeg) )
        cmap = brewmaps('greens', 4);
        cmap = cmap(1:3,:);

        plot(datadiff(indexPos), cccPos(indexPos), 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        x = log10(datadiff(indexPos));
        y = log10(cccPos(indexPos));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        plot(datadiff(indexNeg), cccNeg(indexNeg), 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        x = log10(datadiff(indexNeg));
        y = log10(cccPos(indexNeg));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Sup: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));
    end

    if ( sum(indexNeg) )
        legend('EP', 'SP');
    else
        legend('EP');
    end

    set(gca,'xtick', xtickDiff{i}, 'xticklabel', xtickDiff{i});
    set(gca,'ytick', ytickCCC{i}, 'yticklabel', ytickCCC{i});
    tickpref;
    box off;

    xlim(xaxCCC{i});
    ylim(yaxCCC{i});

    xlabel(xlabelCCC{i});
    ylabel('Correlation Coefficient');
    set(gca,'yscale', 'log');
    set(gca,'xscale', 'log');
    subplot_label(gca,'C');

    set(gcf,'position', [1224 262 321 693]);
    print_mfilename(mfilename);



   fprintf('\n');
   fprintf('%s Larger vs. Smaller\n', datatype{i});

   [r,p] = corrcoef( log10(larger), log10(smaller) );

   fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
   fprintf('%s diff\n', datatype{i});
   x = datadiff(indexPos);

   fprintf('N=%.0f\nMN=%.4f\nSD=%.4f\nSE=%.4f\n', ...
      length(x), mean(x), std(x), std(x)./sqrt(length(x)));
   fprintf('MD=%.4f\nMAD=%.4f\nMX=%.4f\n', ...
      median(datadiff), mad(datadiff), max(datadiff));
    fprintf('\n\n');

end % for i

return;




function iccp_plot_bmf_ccc_exc_excsup_peak(data)


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);

ccpairs_btmf = data.btmf;
ccpairs_bsmf = data.bsmf;


if ( sum(indexNeg) )

    spktype_label{1} = 'CCPairs';
    spktype_xlabel{1} = 'Neuron 1';
    spktype_ylabel{1} = 'Neuron 2';
    datatype = {'btmf', 'bsmf'};

    xaxScatter =  {[12.5 400],[0.05 2.1]};
    yaxScatter = xaxScatter;

    xscale = {'log', 'log', 'log'};
    yscale = {'log', 'log', 'log'};

    edges = {linspace(0, 300, 21), linspace(0, 2, 21)};

    xtickScatter = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
    xtickHist = {[0:60:300], [0:0.4:2]}; 
    xtickDiff = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1 2]}; 


    xlabelScatter = {'bTMF (Hz)', 'bSMF (cyc/oct)'};
    xlabelHist = {'bTMF Difference (Hz)', 'bSMF Difference (cyc/oct)'};
    xlabelCCC = xlabelHist;

    yaxDiff = {[0 0.4], [0 0.6]};

    ytickHist = {0:0.1:0.4, 0:0.2:0.6, linspace(0,0.6,4)}; 

    ytickCCC = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1]}; 


    xaxCCC = {[0.01 400],[0.001 2]};
    yaxCCC = {[0.001 1], [0.001 1]};

    markersize = 4;



    for i = 1:length(datatype)


        figure;
        dataspktype = eval(['ccpairs_' datatype{i}]);
        dataspktype_pos = dataspktype(indexPos,:);
        dataspktype_pos_neg = dataspktype(indexPosNeg,:);

        [larger,smaller] = iccp_largersmaller(dataspktype(:,1),dataspktype(:,2));
        datadiff = abs( larger - smaller );


        [larger_pos,smaller_pos] = ...
            iccp_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));
        datadiff_pos = abs( larger_pos - smaller_pos );


        [larger_pos_neg,smaller_pos_neg] = ...
            iccp_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));
        datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );



        % Make a random distribution, sample, compute paired differences,
        % and plot.
        dataUnq = unique(dataspktype);
        indexRand = ceil(length(dataUnq) * rand(1,size(dataspktype,1)));
        sampleRand1 = dataUnq(indexRand);


        indexRand = ceil(length(dataUnq) * rand(1,size(dataspktype,1)));
        sampleRand2 = dataUnq(indexRand);

        [larger,smaller] = iccp_largersmaller(sampleRand1,sampleRand2);
        datadiffRand = abs( larger - smaller );


pval = ranksum(datadiffRand, datadiff_pos)
pval = ranksum(datadiffRand, datadiff_pos_neg)



        subplot(3,1,1);
        hold on;
        plot([xaxScatter{i}], [yaxScatter{i}], 'k-');


        cmap = brewmaps('blues', 4);
        cmap = cmap(1:3,:);

        plot(larger_pos, smaller_pos, 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        plot(larger_pos_neg, smaller_pos_neg, 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
        tickpref;
        set(gca,'xtick', xtickScatter{i}, 'xticklabel', xtickScatter{i});
        set(gca,'ytick', xtickScatter{i}, 'yticklabel', xtickScatter{i});

        xlim([xaxScatter{i}]);
        ylim([yaxScatter{i}]);
        xlabel(sprintf('%s %s', xlabelScatter{i}, spktype_xlabel{1}));
        ylabel(sprintf('%s %s', xlabelScatter{i}, spktype_ylabel{1}));
        subplot_label(gca,'A');




        subplot(3,1,2);
        hold on;

        n = histc(datadiff_pos, edges{i});
        pdf_pos = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos, 's-', 'markersize', markersize, ...
            'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
        set(hp, 'color', cmap(1,:));

        n = histc(datadiff_pos_neg, edges{i});
        pdf_pos_neg = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos_neg, 's-', 'markersize', markersize, ...
            'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
        set(hp, 'color', cmap(3,:));


        n = histc(datadiffRand, edges{i});
        pdfRand = n ./ sum(n);
        hp = plot(edges{i}, pdfRand, 'k-');
        set(hp, 'color', 'k', 'linewidth', 2);


%         hp = plot(edges{i}, pdfRand, 's-', ...
%             'markersize', markersize, ...
%             'markerfacecolor', 'k', ...
%             'markeredgecolor', 'k' );
%         set(hp, 'color', 'k');


        hold on;
        box off;
        tickpref;
        set(gca,'xtick', xtickHist{i}, 'xticklabel', xtickHist{i});
        set(gca,'ytick', ytickHist{i}, 'yticklabel', ytickHist{i});
        range = max(edges{i}) - min(edges{i});
        xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
        ylim([yaxDiff{i}]);

        xlabel(xlabelHist{i});
        ylabel('Proportion');

        legend('EP', 'EP+SP', 'Rand');
        subplot_label(gca,'B');



        subplot(3,1,3);
        hold on;

        cmap = brewmaps('greens', 4);
        cmap = cmap(1:3,:);

        plot(datadiff(indexPos), cccPos(indexPos), 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        x = log10(datadiff(indexPos));
        y = log10(cccPos(indexPos));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        plot(datadiff(indexNeg), cccNeg(indexNeg), 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        x = log10(datadiff(indexNeg));
        y = log10(cccPos(indexNeg));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Sup: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        legend('EP', 'SP');

        set(gca,'xtick', xtickDiff{i}, 'xticklabel', xtickDiff{i});
        set(gca,'ytick', ytickCCC{i}, 'yticklabel', ytickCCC{i});
        tickpref;
        box off;

        xlim(xaxCCC{i});
        ylim(yaxCCC{i});

        xlabel(xlabelCCC{i});
        ylabel('Correlation Coefficient');
        set(gca,'yscale', 'log');
        set(gca,'xscale', 'log');
        subplot_label(gca,'C');

        set(gcf,'position', [1224 262 321 693]);
        print_mfilename(mfilename);



        fprintf('\n');
        fprintf('%s Larger vs. Smaller\n', datatype{i});

        [r,p] = corrcoef( log10(larger), log10(smaller) );

        fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
        fprintf('%s diff\n', datatype{i});
        x = datadiff(indexPos);

        fprintf('N=%.0f\nMN=%.4f\nSD=%.4f\nSE=%.4f\n', ...
            length(x), mean(x), std(x), std(x)./sqrt(length(x)));
        fprintf('MD=%.4f\nMAD=%.4f\nMX=%.4f\n', ...
            median(datadiff), mad(datadiff), max(datadiff));
        fprintf('\n\n');

    end % (for i)

end % (if)

return;




function iccp_plot_bmf_ccc_exc_excsup_sup_peak(data)

cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);

ccpairs_btmf = data.btmf;
ccpairs_bsmf = data.bsmf;


if ( sum(indexNeg) )

    spktype_label{1} = 'CCPairs';
    spktype_xlabel{1} = 'Neuron 1';
    spktype_ylabel{1} = 'Neuron 2';
    datatype = {'btmf', 'bsmf'};

    xaxScatter =  {[12.5 400],[0.05 2.1]};
    yaxScatter = xaxScatter;

    xscale = {'log', 'log', 'log'};
    yscale = {'log', 'log', 'log'};

    edges = {linspace(0, 300, 21), linspace(0, 2, 21)};

    xtickScatter = {[12.5 25 50 100 200 400], [0.01 0.1 1 2]};
    xtickHist = {[0:60:300], [0:0.4:2]}; 
    xtickDiff = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1 2]}; 


    xlabelScatter = {'bTMF (Hz)', 'bSMF (cyc/oct)'};
    xlabelHist = {'bTMF Difference (Hz)', 'bSMF Difference (cyc/oct)'};
    xlabelCCC = xlabelHist;

    yaxDiff = {[0 0.4], [0 0.6]};

    ytickHist = {0:0.1:0.4, 0:0.2:0.6, linspace(0,0.6,4)}; 

    ytickCCC = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1]}; 


    xaxCCC = {[0.01 400],[0.001 2]};
    yaxCCC = {[0.001 1], [0.001 1]};

    markersize = 4;


    for i = 1:length(datatype)

        dataspktype = eval(['ccpairs_' datatype{i}]);
        dataspktype_pos = dataspktype(indexPos,:);
        dataspktype_neg = dataspktype(indexNeg,:);
        dataspktype_pos_neg = dataspktype(indexPosNeg,:);

        [larger,smaller] = iccp_largersmaller(dataspktype(:,1),dataspktype(:,2));
        datadiff = abs( larger - smaller );

        [larger_pos,smaller_pos] = ...
            iccp_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));
        datadiff_pos = abs( larger_pos - smaller_pos );

        [larger_neg,smaller_neg] = ...
            iccp_largersmaller(dataspktype_neg(:,1),dataspktype_neg(:,2));
        datadiff_neg = abs( larger_neg - smaller_neg );

        [larger_pos_neg,smaller_pos_neg] = ...
            iccp_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));
        datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );
        iccp_ccc_pos_neg_summary(datadiff_pos, datadiff_neg, datadiff_pos_neg, datatype{i});


        figure;

        subplot(3,1,1);
        hold on;
        plot([xaxScatter{i}], [yaxScatter{i}], 'k-');

        cmap = brewmaps('blues', 4);
        cmap = cmap(1:3,:);

        plot(larger_pos, smaller_pos, 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        plot(larger_neg, smaller_neg, 's', 'color', cmap(2,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(2,:), ...
        'markeredgecolor', cmap(2,:));

        plot(larger_pos_neg, smaller_pos_neg, 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
        tickpref;
        set(gca,'xtick', xtickScatter{i}, 'xticklabel', xtickScatter{i});
        set(gca,'ytick', xtickScatter{i}, 'yticklabel', xtickScatter{i});

        xlim([xaxScatter{i}]);
        ylim([yaxScatter{i}]);
        xlabel(sprintf('%s %s', xlabelScatter{i}, spktype_xlabel{1}));
        ylabel(sprintf('%s %s', xlabelScatter{i}, spktype_ylabel{1}));
        subplot_label(gca,'A');



        subplot(3,1,2);
        hold on;

        n = histc(datadiff_pos, edges{i});
        pdf_pos = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos, 's-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
        set(hp, 'color', cmap(1,:));

        n = histc(datadiff_neg, edges{i});
        pdf_neg = n ./ sum(n);
        hp = plot(edges{i}, pdf_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(2,:), 'markeredgecolor', cmap(2,:) );
        set(hp, 'color', cmap(2,:));

        n = histc(datadiff_pos_neg, edges{i});
        pdf_pos_neg = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
        set(hp, 'color', cmap(3,:));

        box off;
        tickpref;
        set(gca,'xtick', xtickHist{i}, 'xticklabel', xtickHist{i});
        set(gca,'ytick', ytickHist{i}, 'yticklabel', ytickHist{i});
        range = max(edges{i}) - min(edges{i});
        xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
        ylim([yaxDiff{i}]);

        xlabel(xlabelHist{i});
        ylabel('Proportion');

        legend('EP', 'SP','EP+SP');
        subplot_label(gca,'B');



        subplot(3,1,3);
        hold on;

        cmap = brewmaps('greens', 4);
        cmap = cmap(1:3,:);

        plot(datadiff(indexPos), cccPos(indexPos), 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        x = log10(datadiff(indexPos));
        y = log10(cccPos(indexPos));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        plot(datadiff(indexNeg), cccNeg(indexNeg), 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        x = log10(datadiff(indexNeg));
        y = log10(cccPos(indexNeg));
        index = ~isinf(x);
        x = x(index);
        y = y(index);
        [r,p] = corrcoef( x, y );
        fprintf('%s: Diff vs. CCC: Sup: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));
        legend('EP', 'SP');

        set(gca,'xtick', xtickDiff{i}, 'xticklabel', xtickDiff{i});
        set(gca,'ytick', ytickCCC{i}, 'yticklabel', ytickCCC{i});
        tickpref;
        box off;

        xlim(xaxCCC{i});
        ylim(yaxCCC{i});

        xlabel(xlabelCCC{i});
        ylabel('Correlation Coefficient');
        set(gca,'yscale', 'log');
        set(gca,'xscale', 'log');
        subplot_label(gca,'C');

        set(gcf,'position', [1224 262 321 693]);
        print_mfilename(mfilename);



        fprintf('\n');
        fprintf('%s Larger vs. Smaller\n', datatype{i});

        [r,p] = corrcoef( log10(larger), log10(smaller) );

        fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
        fprintf('%s diff\n', datatype{i});
        x = datadiff(indexPos);

        fprintf('N=%.0f\nMN=%.4f\nSD=%.4f\nSE=%.4f\n', ...
        length(x), mean(x), std(x), std(x)./sqrt(length(x)));
        fprintf('MD=%.4f\nMAD=%.4f\nMX=%.4f\n', ...
        median(datadiff), mad(datadiff), max(datadiff));
        fprintf('\n\n');

    end % for i

end % (if)

return;
















