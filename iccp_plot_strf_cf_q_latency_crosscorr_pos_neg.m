function varargout = iccp_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs) 
% iccp_plot_strf_cf_q_latency_crosscorr_pos_neg STRF similarity and correlation strength
% 
%    iccp_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs)
% 
%    ccpairs : struct array. Each element is one ICC pair of neurons. Data for
%    STRFs and the correlation parameters for the pair are stored in the
%    fields of each element of the struct array.
% 
%    The function plots cf, q, and latency for pairs of ICC neurons. These
%    parameters were derived from the STRF. The first plot for either variable 
%    is a scatter plot of the value for one neuron in a pair against the 
%    value for the other neuron. The second plot is a histogram of the 
%    differences between the valus for the pairs of neurons. The third plot 
%    compares the cross-correlation coefficient to the parameter differences 
%    shown in the second plot.
%
%    Two plots are made for each parameter. One is the scatter plot for the
%    parameter. Each point is one pair of neurons, with values for each 
%    neuron on the different axes. The second plot is the difference of the
%    parameter values for each pair of neurons.
%
%    data = iccp_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs) returns
%    the data that were used to make the plots. data is a struct with the
%    following fields: 'sigPos', 'sigNeg', 'pdPos', 'hwPos', 'cccPos', 
%    'pdNeg', 'hwNeg', 'cccNeg', 'cf', 'q', 'latency'.
%
%    pd*, hw*, and ccc* are the peak delay, halfwidth, cross-correlation 
%    coefficient. cf, q, and latency are pure tone parameters. They are 
%    Nx2 arrays, with the data for one pair on each row. sigPos and sigNeg 
%    are vectors of 0's and 1's specifying if the covariance function had 
%    a significant feature for that pair of neurons.


library('dataviz'); % used to add jitter to data points


close all;

data = iccp_ccpairs2ptparams(ccpairs);

%iccp_plot_ptparams_exc_peak(data);

iccp_plot_ptparams_exc_excsup_peak(data);

%iccp_plot_ptparams_exc_excsup_sup_peak(data);


if ( nargout == 1 )
    varargout{1} = data;
end


return;






function data = iccp_ccpairs2ptparams(ccpairs)


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




cf1 = [ccpairs.cf1];
cf2 = [ccpairs.cf2];

q1 = [ccpairs.q1];
q2 = [ccpairs.q2];

latency1 = [ccpairs.latency1];
latency2 = [ccpairs.latency2];

cf = [cf1(:) cf2(:)];
q = [q1(:) q2(:)];
latency = [latency1(:) latency2(:)];


if ( min(min(cf)) > 40 ) % It's in Hz - we want it in kHz
   cf = cf ./ 1000;
end


% Paired pure tone parameter data
data.cf = cf;
data.q = q;
data.latency = latency;

return;




function iccp_plot_ptparams_exc_peak(data)


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);


ccpairs_cf = data.cf;
ccpairs_q = data.q;
ccpairs_latency = data.latency;



datatype = {'cf', 'q', 'latency'};

spktype_xlabel{1} = 'Neuron 1';
spktype_ylabel{1} = 'Neuron 2';


xtickScatter = {[0.5 1 2 4 8 16 32], [0:5], [0:4:20]};
ytickScatter = xtickScatter;

xaxScatter = {[0.5 32], [0 4.5], [0 20]};
yaxScatter = xaxScatter;

xscaleScatter = {'log', 'linear', 'linear'};
yscaleScatter = {'log', 'linear', 'linear'};

xlabelScatter = {'BF (kHz)', 'Q', 'Latency (ms)'};



edges{1} = linspace(0, 2, 11);
edges{2} = linspace(0, 4, 11);
edges{3} = linspace(0, 10, 11);
xtickHist = {edges{1}(1:2:end), edges{2}(1:2:end), edges{3}(1:2:end)} ;
ytickHist = {linspace(0,0.8,5), linspace(0,0.5,5), linspace(0,0.8,5)}; 
yaxHist = {[0 0.8], [0 0.5],[0 0.8]};

xlabelHist = {'BF Difference (oct)', 'Q Difference (oct)', 'Latency Difference (ms)'};

label_left = {'A', 'C', 'E'};
label_right = {'B', 'D', 'F'};

markersize = 2;

figure;

for i = 1:length(datatype)

    dataspktype = eval(['ccpairs_' datatype{i}]);

    if ( i > 2 )
        dataspktype(:,1) = jitter( dataspktype(:,1) );
        dataspktype(:,2) = jitter( dataspktype(:,2) );
    end 

    dataspktype_pos = dataspktype(indexPos,:);

    [larger,smaller] = iccp_largersmaller(dataspktype(:,1),dataspktype(:,2));
    [larger_pos,smaller_pos] = iccp_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));

   if ( strcmp(datatype{i}, 'cf') )
      datadiff = abs( log2( larger ./ smaller ) );
      datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
   elseif ( strcmp(datatype{i}, 'latency') )
      datadiff = abs( larger - smaller );
      datadiff_pos = abs( larger_pos - smaller_pos );
   else % it must be 'q'
      datadiff = abs( log2( larger ./ smaller ) );
      datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
   end



    subplot(3,2,(i-1)*2+1);
    markersize = 2;
    hold on;
    plot([xaxScatter{i}], [yaxScatter{i}], 'k-');

    plot(larger, smaller, 's', ...
    'markersize', markersize, ... %'markerfacecolor', cmap(1,:), ...
    'markeredgecolor', 0.25*ones(1,3)); %cmap(1,:));

    set(gca,'xscale', xscaleScatter{i}, 'yscale', yscaleScatter{i});
    tickpref;

    set(gca,'xtick', xtickScatter{i}, 'xticklabel', xtickScatter{i});
    set(gca,'ytick', ytickScatter{i}, 'yticklabel', ytickScatter{i});
    xlim([xaxScatter{i}]);
    ylim([yaxScatter{i}]);

    xlabel(sprintf('%s %s', xlabelScatter{i}, spktype_xlabel{1}));
    ylabel(sprintf('%s %s', xlabelScatter{i}, spktype_ylabel{1}));

    subplot_label(gca, label_left{i});



    subplot(3,2,2*i);
    markersize = 4;
    hold on;

    n = histc(datadiff, edges{i});
    pdf_pos = n ./ sum(n);
    plot(edges{i}, pdf_pos, 'ks-', 'markersize', 5, ...%markersize, ...
        'markerfacecolor', 0.7*ones(1,3), ...
        'markeredgecolor', zeros(1,3)); 

    box off;
    tickpref;


    set(gca,'xtick', xtickHist{i}, 'xticklabel', xtickHist{i});
    set(gca,'ytick', ytickHist{i}, 'yticklabel', ytickHist{i});
    range = max(edges{i}) - min(edges{i});
    xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
    ylim([yaxHist{i}]);

    xlabel(xlabelHist{i});
    ylabel('Proportion');

    legend('EP');
    subplot_label(gca,label_right{i});

end % for i

set(gcf,'position', [724 300 438 583]);
print_mfilename(mfilename);

return;




function iccp_plot_ptparams_exc_excsup_peak(data)


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);


ccpairs_cf = data.cf;
ccpairs_q = data.q;
ccpairs_latency = data.latency;


if ( sum(indexNeg) )

    datatype = {'cf', 'q', 'latency'};

    spktype_xlabel{1} = 'Neuron 1';
    spktype_ylabel{1} = 'Neuron 2';


    xtickScatter = {[0.5 1 2 4 8 16 32], [0:5], [0:4:20]};
    ytickScatter = xtickScatter;

    xaxScatter = {[0.5 32], [0 4.5], [0 20]};
    yaxScatter = xaxScatter;

    xscaleScatter = {'log', 'linear', 'linear'};
    yscaleScatter = {'log', 'linear', 'linear'};

    xlabelScatter = {'BF (kHz)', 'Q', 'Latency (ms)'};


    edges{1} = linspace(0, 2, 11);
    edges{2} = linspace(0, 4, 11);
    edges{3} = linspace(0, 10, 11);
    xtickHist = {edges{1}(1:2:end), edges{2}(1:2:end), edges{3}(1:2:end)} ;
    ytickHist = {linspace(0,0.8,5), linspace(0,0.5,5), linspace(0,0.8,5)}; 
    yaxHist = {[0 0.8], [0 0.5],[0 0.8]};

    xlabelHist = {'BF Difference (oct)', 'Q Difference (oct)', 'Latency Difference (ms)'};

    label_left = {'A', 'C', 'E'};
    label_right = {'B', 'D', 'F'};

    figure;

    for i = 1:length(datatype)

        dataspktype = eval(['ccpairs_' datatype{i}]);

        if ( i > 2 )
            dataspktype(:,1) = jitter( dataspktype(:,1) );
            dataspktype(:,2) = jitter( dataspktype(:,2) );
        end 

        dataspktype_pos = dataspktype(indexPos,:);
        dataspktype_pos_neg = dataspktype(indexPosNeg,:);

        [larger,smaller] = iccp_largersmaller(dataspktype(:,1),dataspktype(:,2));

        [larger_pos,smaller_pos] = ...
            iccp_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));

        [larger_pos_neg,smaller_pos_neg] = ...
            iccp_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));


        % Make a random distribution, sample, compute paired differences,
        % and plot.
        dataUnq = unique(dataspktype);
        indexRand = ceil(length(dataUnq) * rand(1,size(dataspktype,1)));
        sampleRand1 = dataUnq(indexRand);
        indexRand = ceil(length(dataUnq) * rand(1,size(dataspktype,1)));
        sampleRand2 = dataUnq(indexRand);
        [largerRand,smallerRand] = iccp_largersmaller(sampleRand1,sampleRand2);


        [absc_pos, ord_pos] = ...
            iccp_randomize_columns(dataspktype_pos(:,1),dataspktype_pos(:,2));

        [absc_pos_neg, ord_pos_neg] = ...
            iccp_randomize_columns(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));


        if ( strcmp(datatype{i}, 'cf') )
            datadiff = abs( log2( larger ./ smaller ) );
            datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
            datadiff_pos_neg = abs( log2( larger_pos_neg ./ smaller_pos_neg ) );
            datadiffRand = abs( log2( largerRand ./ smallerRand ) );
        elseif ( strcmp(datatype{i}, 'latency') )
            datadiff = abs( larger - smaller );
            datadiff_pos = abs( larger_pos - smaller_pos );
            datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );

            datadiffRand = abs( larger - smaller );
        else % it must be 'q'
            datadiff = abs( log2( larger ./ smaller ) );
            datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
            datadiff_pos_neg = abs( log2( larger_pos_neg ./ smaller_pos_neg ) );
            datadiffRand = abs( log2( largerRand ./ smallerRand ) );
        end



        logtransform = 1;
        [rpop_med, rpop_ci, pval] = ...
            iccp_pairwise_corr_rand_test(absc_pos, ord_pos, logtransform);
        fprintf('\n');
        fprintf(sprintf('%s Pos pairs - randomized\n', xlabelScatter{i}));
        fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
        fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
        fprintf('\n\n');


        logtransform = 1;
        [rpop_med, rpop_ci, pval] = ...
            iccp_pairwise_corr_rand_test(absc_pos_neg, ord_pos_neg, logtransform);
        fprintf('\n');
        fprintf(sprintf('%s Pos-Neg pairs - randomized\n',xlabelScatter{i}));
        fprintf('Pairwise Randomization: r=%.4f,p=%.4f\n', rpop_med,pval);
        fprintf('Confidence Intervals: [%.4f, %.4f]\n', rpop_ci(1), rpop_ci(2));
        fprintf('\n\n');


pval = ranksum(datadiffRand, datadiff_pos)
pval = ranksum(datadiffRand, datadiff_pos_neg)





        subplot(3,2,(i-1)*2+1);
        markersize = 3;
        hold on;
        plot([xaxScatter{i}], [yaxScatter{i}], 'k-');

        cmap = brewmaps('blues', 4);
        cmap = cmap(1:3,:);

%         plot(larger_pos, smaller_pos, 's', 'color', cmap(1,:), ...
%             'markersize', markersize, 'markerfacecolor', cmap(1,:), ...
%             'markeredgecolor', cmap(1,:));
% 
%         plot(larger_pos_neg, smaller_pos_neg, 's', 'color', cmap(3,:), ...
%             'markersize', markersize, 'markerfacecolor', cmap(3,:), ...
%             'markeredgecolor', cmap(3,:));


        plot(absc_pos, ord_pos, 'o', ...
            'color', 'k', ...
            'markersize', markersize, ...
            'markeredgecolor', 'k');

        gray_face = 0.6 * ones(1,3);
        gray_edge = 0.3 * ones(1,3);

        plot(absc_pos_neg, ord_pos_neg, 's', ...
            'color', cmap(3,:), ...
            'markersize', markersize, ...
            'markerfacecolor', gray_face, ...
            'markeredgecolor', gray_edge);





        set(gca,'xscale', xscaleScatter{i}, 'yscale', yscaleScatter{i});
        tickpref;

        set(gca,'xtick', xtickScatter{i}, 'xticklabel', xtickScatter{i});
        set(gca,'ytick', ytickScatter{i}, 'yticklabel', ytickScatter{i});
        xlim([xaxScatter{i}]);
        ylim([yaxScatter{i}]);

        xlabel(sprintf('%s %s', xlabelScatter{i}, spktype_xlabel{1}));
        ylabel(sprintf('%s %s', xlabelScatter{i}, spktype_ylabel{1}));

        subplot_label(gca, label_left{i});



        subplot(3,2,2*i);
        markersize = 4;
        hold on;

        n = histc(datadiff_pos, edges{i});
        pdf_pos = n ./ sum(n);

%         hp = plot(edges{i}, pdf_pos, 's-', 'markersize', markersize, ...
%             'markerfacecolor', cmap(1,:), 'markeredgecolor', cmap(1,:) );
%         set(hp, 'color', cmap(1,:));
% 
%         n = histc(datadiff_pos_neg, edges{i});
%         pdf_pos_neg = n ./ sum(n);
%         hp = plot(edges{i}, pdf_pos_neg, 's-', 'markersize', markersize, ...
%             'markerfacecolor', cmap(3,:), 'markeredgecolor', cmap(3,:) );
%         set(hp, 'color', cmap(3,:));


        hp = plot(edges{i}, pdf_pos, 'ko-', ...
            'markersize', markersize, ...
            'markerfacecolor', 'k');
        set(hp, 'color', 'k');


        n = histc(datadiff_pos_neg, edges{i});
        pdf_pos_neg = n ./ sum(n);

        hp = plot(edges{i}, pdf_pos_neg, 's-', ...
            'markersize', markersize, ...
            'markerfacecolor', 0.6*ones(1,3), ...
            'markeredgecolor', 0.3*ones(1,3));
        set(hp, 'color', 0.6*ones(1,3));




        % Random distribution
        n = histc(datadiffRand, edges{i});
        pdfRand = n ./ sum(n);
        hp = plot(edges{i}, pdfRand, 'k-');
        set(hp, 'color', 'k', 'linewidth', 1);


        box off;
        tickpref;

        set(gca,'xtick', xtickHist{i}, 'xticklabel', xtickHist{i});
        set(gca,'ytick', ytickHist{i}, 'yticklabel', ytickHist{i});
        range = max(edges{i}) - min(edges{i});
        xlim([min(edges{i})-0.025*range max(edges{i})+0.025*range]);
        ylim([yaxHist{i}]);

        xlabel(xlabelHist{i});
        ylabel('Proportion');

        legend('EP', 'EP+SP');
        subplot_label(gca,label_right{i});

    end % for i

    set(gcf,'position', [724 300 438 583]);
    print_mfilename(mfilename);

end % (if)

return;




function iccp_plot_ptparams_exc_excsup_sup_peak(data)


cccPos = data.cccPos;
cccNeg = data.cccNeg;

indexPos = data.sigPos;
indexNeg = data.sigNeg;
indexPosNeg = indexPos & indexNeg;

cccPosSig = cccPos(indexPos);
cccNegSig = cccNeg(indexNeg);


ccpairs_cf = data.cf;
ccpairs_q = data.q;
ccpairs_latency = data.latency;


if ( sum(indexNeg) )


    datatype = {'cf', 'q', 'latency'};

    spktype_xlabel{1} = 'Neuron 1';
    spktype_ylabel{1} = 'Neuron 2';


    xtickScatter = {[0.5 1 2 4 8 16 32], [0:5], [0:4:20]};
    ytickScatter = xtickScatter;

    xaxScatter = {[0.5 32], [0 4.5], [0 20]};
    yaxScatter = xaxScatter;

    xscaleScatter = {'log', 'linear', 'linear'};
    yscaleScatter = {'log', 'linear', 'linear'};

    xlabelScatter = {'BF (kHz)', 'Q', 'Latency (ms)'};


    edges{1} = linspace(0, 2, 11);
    edges{2} = linspace(0, 4, 11);
    edges{3} = linspace(0, 10, 11);
    xtickHist = {edges{1}(1:2:end), edges{2}(1:2:end), edges{3}(1:2:end)} ;
    ytickHist = {linspace(0,0.8,5), linspace(0,0.5,5), linspace(0,0.8,5)}; 
    yaxHist = {[0 0.8], [0 0.5],[0 0.8]};

    xlabelHist = {'BF Difference (oct)', 'Q Difference (oct)', 'Latency Difference (ms)'};

    label_left = {'A', 'C', 'E'};
    label_right = {'B', 'D', 'F'};

    figure;

    for i = 1:length(datatype)

        dataspktype = eval(['ccpairs_' datatype{i}]);

        if ( i > 2 )
            dataspktype(:,1) = jitter( dataspktype(:,1) );
            dataspktype(:,2) = jitter( dataspktype(:,2) );
        end

        dataspktype_pos = dataspktype(indexPos,:);
        dataspktype_neg = dataspktype(indexNeg,:);
        dataspktype_pos_neg = dataspktype(indexPosNeg,:);

        [larger,smaller] = iccp_largersmaller(dataspktype(:,1),dataspktype(:,2));

        [larger_pos,smaller_pos] = iccp_largersmaller(dataspktype_pos(:,1),dataspktype_pos(:,2));

        [larger_neg,smaller_neg] = ...
            iccp_largersmaller(dataspktype_neg(:,1),dataspktype_neg(:,2));

        [larger_pos_neg,smaller_pos_neg] = ...
            iccp_largersmaller(dataspktype_pos_neg(:,1),dataspktype_pos_neg(:,2));


        if ( strcmp(datatype{i}, 'cf') )
            datadiff = abs( log2( larger ./ smaller ) );
            datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
            datadiff_neg = abs( log2( larger_neg ./ smaller_neg ) );
            datadiff_pos_neg = abs( log2( larger_pos_neg ./ smaller_pos_neg ) );
        elseif ( strcmp(datatype{i}, 'latency') )
            datadiff = abs( larger - smaller );
            datadiff_pos = abs( larger_pos - smaller_pos );
            datadiff_neg = abs( larger_neg - smaller_neg );
            datadiff_pos_neg = abs( larger_pos_neg - smaller_pos_neg );
        else % it must be 'q'
            datadiff = abs( log2( larger ./ smaller ) );
            datadiff_pos = abs( log2( larger_pos ./ smaller_pos ) );
            datadiff_neg = abs( log2( larger_neg ./ smaller_neg ) );
            datadiff_pos_neg = abs( log2( larger_pos_neg ./ smaller_pos_neg ) );
        end

        iccp_ccc_pos_neg_summary(datadiff_pos, datadiff_neg, datadiff_pos_neg, datatype{i});



        subplot(3,2,(i-1)*2+1);
        markersize = 2;
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

        set(gca,'xscale', xscaleScatter{i}, 'yscale', yscaleScatter{i});

        tickpref;

        set(gca,'xtick', xtickScatter{i}, 'xticklabel', xtickScatter{i});
        set(gca,'ytick', ytickScatter{i}, 'yticklabel', ytickScatter{i});
        xlim([xaxScatter{i}]);
        ylim([yaxScatter{i}]);

        xlabel(sprintf('%s %s', xlabelScatter{i}, spktype_xlabel{1}));
        ylabel(sprintf('%s %s', xlabelScatter{i}, spktype_ylabel{1}));

        subplot_label(gca, label_left{i});



        subplot(3,2,2*i);
        markersize = 4;
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
        ylim([yaxHist{i}]);

        xlabel(xlabelHist{i});
        ylabel('Proportion');

        legend('EP', 'SP','EP+SP');
        subplot_label(gca,label_right{i});

    end % for i

    set(gcf,'position', [724 300 438 583]);
    print_mfilename(mfilename);

end % (if)

return;















