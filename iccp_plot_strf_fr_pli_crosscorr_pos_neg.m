function varargout = iccp_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs)
% iccp_plot_strf_fr_pli_crosscorr_pos_neg STRF similarity and correlation strength
% 
%    iccp_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs)
% 
%    ccpairs : struct array. Each element is one ICC pair of neurons. Data for
%    STRFs and the correlation parameters for the pair are stored in the
%    fields of each element of the struct array.
% 
%    The function plots firing rate and phase-locking index values for pairs
%    of ICC neurons. The first plot for either variable is a scatter plot of 
%    the value for one neuron in a pair against the value for the other neuron.
%    The second plot is a histogram of the differences between the valus for
%    the pairs of neurons. The third plot compares the cross-correlation
%    coefficient to the parameter differences shown in the second plot.
%
%    These three plots are made for each parameter.
%
%    Data is broken into multiple groups. The CCC data is broken into two
%    groups: the strength of the positive cross-covariance function peak, and
%    the negative peak. Both strengths are plotted against STRF similarity.
%    Only values for significant peaks are plotted.
%
%    data = iccp_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs) returns
%    the data that were used to make the plots. data is a struct with the
%    following fields: 'sigPos', 'sigNeg', 'pdPos', 'hwPos', 'cccPos', 
%    'pdNeg', 'hwNeg', 'cccNeg', 'fr', 'pli', 'spi'.
%
%    pd*, hw*, and ccc* are the peak delay, halfwidth, cross-correlation 
%    coefficient. fr, pli, and spi are the firing rate, phase-locking index,
%    and the strf separability. They are Nx2 arrays, with the data for one
%    pair on each row. sigPos and sigNeg are vectors of 0's and 1's 
%    specifying if the covariance function had a significant feature for 
%    that pair of neurons.


close all;

data = iccp_get_ccc_fr_pli(ccpairs);

% iccp_plot_ccc_fr_pli_exc_peak(data);

iccp_plot_ccc_fr_pli_exc_excsup_peak(data);

%iccp_plot_ccc_fr_pli_exc_excsup_sup_peak(data);


if ( nargout == 1 )
    varargout{1} = data;
end


return;





function data = iccp_get_ccc_fr_pli(ccpairs)


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


fr1 = [ccpairs.fr1];
fr2 = [ccpairs.fr2];

pli1 = [ccpairs.pli1];
pli2 = [ccpairs.pli2];

spi1 = [ccpairs.sepindex1];
spi2 = [ccpairs.sepindex2];

data.fr = [fr1(:) fr2(:)];
data.pli = [pli1(:) pli2(:)];
data.spi = [spi1(:) spi2(:)];


return;




function iccp_plot_ccc_fr_pli_exc_peak(data)

ccpairs_fr = data.fr;
ccpairs_pli = data.pli;
ccpairs_si = data.spi;


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);




% Set up plotting parameters, labels, etc.
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

ymaxdiff = [0.75 0.6 0.65];
ytick = {0:0.1:0.7, linspace(0,0.6,5), linspace(0,0.6,4)}; 
xtick2 = {0:20:100, 0:0.1:0.4, [0.0001 0.001 0.01 0.1 1]}; 
xtick3 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
yax = {[0.001 1], [0.001 1], [0.0001 2]};
xax = {[0.01 100],[0.0001 1] [0.0001 2]};

markersize = 2;

for i = 1:2 %length(datatype)


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



    subplot(3,1,1);
    hold on;
    plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');


    plot(larger_pos, smaller_pos, 's', ...
    'markersize', markersize, ... %'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', 0.25*ones(1,3)); %cmap(1,:));


    set(gca,'xscale', xscale{i}, 'yscale', yscale{i});
    tickpref;
    set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
    set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
    xlim([xmin(i) xmax(i)]);
    ylim([ymin(i) ymax(i)]);
    xlabel( dataxlabel{i} );
    ylabel( dataylabel{i} );
    subplot_label(gca,'A');



    subplot(3,1,2);
    hold on;

    n = histc(datadiff_pos, edges{i});
    pdf_pos = n ./ sum(n);
    plot(edges{i}, pdf_pos, 'ks-', 'markersize', 5, ...%markersize, ...
        'markerfacecolor', 0.5*ones(1,3), ...
        'markeredgecolor', zeros(1,3)); 

    hold on;
    box off;
    tickpref;
    set(gca,'xtick', xtick2{i}, 'xticklabel', xtick2{i});
    set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
    range = max(edges{i}) - min(edges{i});
    xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
    ylim([0 ymaxdiff(i)]);
    ylabel('Proportion');
    xlabel(xlabel2{i});
    legend('EP');
    subplot_label(gca,'B');


     
    subplot(3,1,3);
    hold on;

    plot(datadiff(indexPos), cccPos(indexPos), 's', ...
    'markersize', markersize, ... %'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', 0.25*ones(1,3)); %cmap(1,:));

    [r,p] = corrcoef( log10(datadiff(indexPos)), log10(cccPos(indexPos)) );
    fprintf('%s: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

    legend('EP');

    set(gca,'xtick', xtick3{i}, 'xticklabel', xtick3{i});
    set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
    tickpref;
    box off;
    xlim(xax{i});
    ylim(yax{i});
    xlabel(xlabel2{i});
    ylabel('Correlation Coefficient');
    set(gca,'yscale', 'log');
    set(gca,'xscale', 'log');
    subplot_label(gca,'C');


    fprintf('\n');
    fprintf('%s\n', datatype{i});
    [r,p] = corrcoef( log10(larger), log10(smaller) );
    fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
    fprintf('%s diff\n', datatype{i});
    x = datadiff(indexPos);
    fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
        length(x), mean(x), std(x), std(x)./sqrt(length(x)));
    fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
        median(datadiff), mad(datadiff), max(datadiff));
    x = datadiff(indexPos);
    y = cccPos(indexPos);
    [r,p] = corrcoef(x(:), log10(y(:)));
    fprintf('ccc vs diff, r = %.4f, p = %.4f', r(2), p(2));
    fprintf('\n');

    set(gcf,'position', [1434 131 321 694]);
    print_mfilename(mfilename);

end % for i


return;




function iccp_plot_ccc_fr_pli_exc_excsup_peak(data)

ccpairs_fr = data.fr;
ccpairs_pli = data.pli;
ccpairs_si = data.spi;

cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);


if ( sum(indexNeg) )


    % Set up plotting parameters, labels, etc.
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

    ymaxdiff = [0.75 0.6 0.65];
    ytick = {0:0.1:0.7, linspace(0,0.6,5), linspace(0,0.6,4)}; 
    xtick2 = {0:20:100, 0:0.1:0.4, [0.0001 0.001 0.01 0.1 1]}; 
    xtick3 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
    ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
    yax = {[0.001 1], [0.001 1], [0.0001 2]};
    xax = {[0.01 100],[0.0001 1] [0.0001 2]};

    markersize = 4;

    for i = 1:2 %length(datatype)


        figure;

        % Get data over population, for positive peaks, and negative peaks
        dataspktype = eval(['ccpairs_' datatype{i}]);
        dataspktype_pos = dataspktype(indexPos,:);
        dataspktype_pos_neg = dataspktype(indexPosNeg,:);


        % Order data and take differences
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
        plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');


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
        set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
        set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
        xlim([xmin(i) xmax(i)]);
        ylim([ymin(i) ymax(i)]);
        xlabel( dataxlabel{i} );
        ylabel( dataylabel{i} );
        subplot_label(gca,'A');



        subplot(3,1,2);
        hold on;

        n = histc(datadiff_pos, edges{i});
        pdf_pos = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos, 's-', 'markersize', markersize, 'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
        set(hp, 'color', cmap(1,:));

        n = histc(datadiff_pos_neg, edges{i});
        pdf_pos_neg = n ./ sum(n);
        hp = plot(edges{i}, pdf_pos_neg, 's-', 'markersize', markersize, 'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
        set(hp, 'color', cmap(3,:));

%         n = histc(datadiffRand, edges{i});
%         pdfRand = n ./ sum(n);
%         hp = plot(edges{i}, pdfRand, 's-', 'markersize', markersize, 'markerfacecolor', 'k', 'markeredgecolor', 'k' );
%         set(hp, 'color', 'k');

        n = histc(datadiffRand, edges{i});
        pdfRand = n ./ sum(n);
        hp = plot(edges{i}, pdfRand, 'k-');
        set(hp, 'color', 'k', 'linewidth', 2);


        hold on;
        box off;
        tickpref;
        set(gca,'xtick', xtick2{i}, 'xticklabel', xtick2{i});
        set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
        range = max(edges{i}) - min(edges{i});
        xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
        ylim([0 ymaxdiff(i)]);
        ylabel('Proportion');
        xlabel(xlabel2{i});
        legend('EP', 'EP+SP', 'Rand');
        subplot_label(gca,'B');


        subplot(3,1,3);
        hold on;

        cmap = brewmaps('greens', 4);
        cmap = cmap(1:3,:);

        plot(datadiff(indexPos), cccPos(indexPos), 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        [r,p] = corrcoef( log10(datadiff(indexPos)), log10(cccPos(indexPos)) );
        fprintf('%s: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        plot(datadiff(indexNeg), cccNeg(indexNeg), 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        [r,p] = corrcoef( log10(datadiff(indexNeg)), log10(cccNeg(indexNeg)) );
        fprintf('%s: Sup: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        legend('EP', 'SP');

       set(gca,'xtick', xtick3{i}, 'xticklabel', xtick3{i});
       set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
       tickpref;
       box off;
       xlim(xax{i});
       ylim(yax{i});
       xlabel(xlabel2{i});
       ylabel('Correlation Coefficient');
       set(gca,'yscale', 'log');
       set(gca,'xscale', 'log');
       subplot_label(gca,'C');


       fprintf('\n');
       fprintf('%s\n', datatype{i});
       [r,p] = corrcoef( log10(larger), log10(smaller) );
       fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
       fprintf('%s diff\n', datatype{i});
       x = datadiff(indexPos);

       fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
          length(x), mean(x), std(x), std(x)./sqrt(length(x)));
       fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
          median(datadiff), mad(datadiff,1), max(datadiff));

       x = datadiff(indexPos);
       y = cccPos(indexPos);
       [r,p] = corrcoef(x(:), log10(y(:)));
       fprintf('ccc vs diff, r = %.4f, p = %.4f', r(2), p(2));
       fprintf('\n');

       set(gcf,'position', [1434 131 321 694]);
       print_mfilename(mfilename);

    end % for i

end % (if)


return;




function iccp_plot_ccc_fr_pli_exc_excsup_sup_peak(data)

ccpairs_fr = data.fr;
ccpairs_pli = data.pli;
ccpairs_si = data.spi;


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);


if ( sum(indexNeg) )

    % Set up plotting parameters, labels, etc.
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

    ymaxdiff = [0.75 0.6 0.65];
    ytick = {0:0.1:0.7, linspace(0,0.6,5), linspace(0,0.6,4)}; 
    xtick2 = {0:20:100, 0:0.1:0.4, [0.0001 0.001 0.01 0.1 1]}; 
    xtick3 = {[0.01 0.1 1 10 100], [0.0001 0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
    ytick2 = {[0.001 0.01 0.1 1], [0.001 0.01 0.1 1], [0.0001 0.001 0.01 0.1 1]}; 
    yax = {[0.001 1], [0.001 1], [0.0001 2]};
    xax = {[0.01 100],[0.0001 1] [0.0001 2]};

    markersize = 2;

    for i = 1:2 %length(datatype)


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

        [larger_neg,smaller_neg] = ...
            iccp_largersmaller(dataspktype_neg(:,1),dataspktype_neg(:,2));
        datadiff_neg = abs( larger_neg - smaller_neg );

        [larger_pos_neg,smaller_pos_neg] = ...
            iccp_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));
        datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );

        iccp_ccc_pos_neg_summary(datadiff_pos, datadiff_neg, datadiff_pos_neg, datatype{i});


        subplot(3,1,1);
        hold on;
        plot([xmin(i) xmax(i)], [ymin(i) ymax(i)], 'k-');

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
        set(gca,'xtick', tick{i}, 'xticklabel', tick{i});
        set(gca,'ytick', tick{i}, 'yticklabel', tick{i});
        xlim([xmin(i) xmax(i)]);
        ylim([ymin(i) ymax(i)]);
        xlabel( dataxlabel{i} );
        ylabel( dataylabel{i} );
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
        set(gca,'xtick', xtick2{i}, 'xticklabel', xtick2{i});
        set(gca,'ytick', ytick{i}, 'yticklabel', ytick{i});
        range = max(edges{i}) - min(edges{i});
        xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
        ylim([0 ymaxdiff(i)]);
        ylabel('Proportion');
        xlabel(xlabel2{i});
        legend('EP', 'SP','EP+SP');
        subplot_label(gca,'B');


        subplot(3,1,3);
        hold on;

        cmap = brewmaps('greens', 4);
        cmap = cmap(1:3,:);

        plot(datadiff(indexPos), cccPos(indexPos), 's', 'color', cmap(1,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
        'markeredgecolor', cmap(1,:));

        [r,p] = corrcoef( log10(datadiff(indexPos)), log10(cccPos(indexPos)) );
        fprintf('%s: Exc: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        plot(datadiff(indexNeg), cccNeg(indexNeg), 's', 'color', cmap(3,:), ...
        'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
        'markeredgecolor', cmap(3,:));

        [r,p] = corrcoef( log10(datadiff(indexNeg)), log10(cccNeg(indexNeg)) );
        fprintf('%s: Sup: r=%.4f, p = %.4f\n', datatype{i}, r(2),p(2));

        legend('EP', 'SP');

        set(gca,'xtick', xtick3{i}, 'xticklabel', xtick3{i});
        set(gca,'ytick', ytick2{i}, 'yticklabel', ytick2{i});
        tickpref;
        box off;
        xlim(xax{i});
        ylim(yax{i});
        xlabel(xlabel2{i});
        ylabel('Correlation Coefficient');
        set(gca,'yscale', 'log');
        set(gca,'xscale', 'log');
        subplot_label(gca,'C');


        fprintf('\n');
        fprintf('%s\n', datatype{i});
        [r,p] = corrcoef( log10(larger), log10(smaller) );
        fprintf('r=%.4f, p = %.4f\n', r(2),p(2));
        fprintf('%s diff\n', datatype{i});
        x = datadiff(indexPos);
        fprintf('N=%.0f, MN=%.4f, SD=%.4f, SE=%.4f\n', ...
            length(x), mean(x), std(x), std(x)./sqrt(length(x)));
        fprintf('MD=%.4f, MAD=%.4f, MX=%.4f\n', ...
            median(datadiff), mad(datadiff), max(datadiff));
        x = datadiff(indexPos);
        y = cccPos(indexPos);
        [r,p] = corrcoef(x(:), log10(y(:)));
        fprintf('ccc vs diff, r = %.4f, p = %.4f', r(2), p(2));
        fprintf('\n');

       set(gcf,'position', [1434 131 321 694]);
       print_mfilename(mfilename);

    end % (for i)

end % (if)

return;

































