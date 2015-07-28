function index = iccpairs_select_synchronous_ccpairs(ccpairs)
%pairs_select_synchronous_ccpairs Identify sharply synchronous neuron pairs
%
%   index = pairs_select_synchronous_ccpairs(ccpairs)
%
%   ccpairs : database of neurons recorded from the same channel. ccpairs
%   is a struct array, with each element holding data for one pair. The
%   fields of ccpairs hold receptive field parameters cross-covariance data
%   for each pair.
%
%   strf : optional argument. It is the strf struct array for the 
%   individual neurons in ccpairs. If strf is provided, then the
%   strfs will be shown for the two constituent neurons in the
%   ccpairs element. 
%
%   trigger : vector of trigger times that is to accompany STRF.
%
%   index : indicates elements of ccpairs that the user has identified as
%   a synchronous spiking pair.
% 
%   index = pairs_select_synchronous_ccpairs(ccpairs)


if ( nargin ~= 1 && nargin ~= 3 )
    error('Wrong input arguments.');
end

if ( exist('temp-ccpairs-select.mat', 'file') )
    load('temp-ccpairs-select.mat')
    istart = max(index) + 1;
else
    index = [];
    istart = 1;
end




figure;

index = [];

for i = istart:length(ccpairs)

    clf;

    exp = ccpairs(i).exp;
    site = ccpairs(i).site;
    chan = ccpairs(i).chan;
    model1 = ccpairs(i).model1;
    model2 = ccpairs(i).model2;
    position = ccpairs(i).position;
    fs = ccpairs(i).fs;

    pd = ccpairs(i).peakdelay;
    hw = ccpairs(i).halfwidth;
    ccc = ccpairs(i).ccc;

    conf_limit = ccpairs(i).conf_limit;
    delay = ccpairs(i).delay;
    q12 = ccpairs(i).q12;
    r12 = ccpairs(i).r12;


    spiketimes1 = ccpairs(i).spiketimes1;
    spiketimes2 = ccpairs(i).spiketimes2;


    nsp1 = length(spiketimes1);
    nsp2 = length(spiketimes2);


%     if ( nargin == 3 )
        % Put spiketimes in seconds
        a = spiketimes1 / 1000;
        b = spiketimes2 / 1000;
        binwidth = 1.5; % bin resolution
        [ccdelay, ccraw, ccflat, ccnorm] = agmon_crosscorr_binless(a, a, binwidth);
%         [ccdelay, ccraw, ccflat, ccnorm] = agmon_crosscorr_binless(a, b, binwidth);
%     end


   % Get spike waveforms
   % ----------------------------------------
   wavetime1 = ccpairs(i).meanwaveform1(:,1);
   waveform1 = ccpairs(i).meanwaveform1(:,2);
   waveform1= waveform1 ./ 2048;

   wavetime2 = ccpairs(i).meanwaveform2(:,1);
   waveform2 = ccpairs(i).meanwaveform2(:,2);
   waveform2= waveform2 ./ 2048;

   % Normalize the waveforms to [-1, 1]
   max1 = max([max(abs(waveform1)) max(abs(min(waveform1)))]);
   max2 = max([max(abs(waveform2)) max(abs(min(waveform2)))]);
   maxmax = max([max1 max2]);
   waveform1 = waveform1 ./ maxmax;
   waveform2 = waveform2 ./ maxmax;

 
   % Plot the waveforms
   % -------------------------------
    if ( nargin == 1 )
        subplot(2,2,1);
    else
        subplot(2,3,1);
    end
    
   hold on;
   plot(wavetime1, waveform1, 'k-', 'linewidth', 1.5);
   plot(wavetime2, waveform2, '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5);
   xlim([-1 2]);
   ylim([-1.05 1.05]);
   box off;
   tickpref;
   set(gca,'xtick', [-1:2], 'xticklabel', [-1:2]);
   xlabel('Time (ms)');
   set(gca,'yticklabel', [-1:1], 'ytick', [-1:1]);
   ylabel('Norm. Amp.');
    title(sprintf('#%.0f of %.0f: %s site%.0f chan%.0f model1[%.0f] model2[%.0f]', ...
    i, length(ccpairs), exp, site, chan, model1, model2));


   % Plot the cross-covariance
   % -------------------------------
    if ( nargin == 1 )
        subplot(2,2,3);
    else
        subplot(2,3,4);
    end
    upper95qab = conf_limit;
    lower95qab = -conf_limit;
    hold on;
    ymin = min([min(q12) lower95qab]);
    ymax = max([max(q12) upper95qab]);
    range = ymax-ymin;
    plot([0 0], [ymin-0.1*range ymax+0.1*range], '-', ...
        'linewidth', 3, 'color', 0.5*ones(1,3));
    plot([min(delay) max(delay)], [0 0], '-', ...
        'linewidth', 3, 'color', 0.5*ones(1,3));
    bar(delay, q12, 'k');

    plot([min(delay) max(delay)], [upper95qab upper95qab], 'r-');
    plot([min(delay) max(delay)], [lower95qab lower95qab], 'r-');

    plot([pd pd], [ymin-0.1*range ymax+0.1*range], 'r-');

    xlim([min(delay) max(delay)]);
    xlim([-25 25]);
    xtick = -20:10:20;
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    xlabel('Delay (ms)');
    ylabel('Cross Cov');
    ylim([ymin-0.1*range ymax+0.1*range]);
    tickpref;
    title(sprintf('#%.0f of %.0f: PD = %.2f, HW = %.2f', ...
        i, length(ccpairs), pd, hw));



    if ( nargin == 1 )
        subplot(2,2,2);
    else
        subplot(2,3,2);
    end
    hold on;
    plot([0 0], [0 1.1*max(r12)], '-', ...
        'linewidth', 3, 'color', 0.5*ones(1,3));
    bar(delay, r12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
    plot([pd pd], [0 1.05*max(r12)], 'r-');
    xlim([min(delay) max(delay)]);
    xlim([-25 25]);
    ylim([0 1.1*max(r12)]);
    xtick = -20:10:20;
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    xlabel('Delay (ms)');
    ylabel('Cross-correlation');
    %    ylim([ymin-0.1*range ymax+0.1*range]);
    tickpref;


    if ( nargin == 3 )
        subplot(2,3,5);
        index_strf = find_strf_element(strf, exp, site, chan, model1);
        plot_strf_single(strf(index_strf),trigger);

        subplot(2,3,6);
        index_strf = find_strf_element(strf, exp, site, chan, model2);
        plot_strf_single(strf(index_strf),trigger);

        cmap = brewmaps('rdbu', 21);
        cmap = cmap([1:8 11 14:end],:);
        colormap(cmap);
    end




    if ( nargin == 1 )
        subplot(2,2,4);
    else ( nargin == 3 )
        subplot(2,3,3);
    end
    hold on;
    plot([0 0], [0 1.1*max(ccnorm)], '-', ...
        'linewidth', 3, 'color', 0.5*ones(1,3));
    hb = bar(ccdelay, ccnorm);
    set(hb,'facecolor', 0*ones(1,3));
    set(hb,'edgecolor', 0*ones(1,3));
    set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
    xlim([-max(ccdelay)-binwidth max(ccdelay)+binwidth]);
    ylim([0 1.1*max(ccnorm)]);
    box off;
    xlabel('Delay (ms)');
    ylabel('CC');
    title(sprintf('Binwidth = %.2f ms', binwidth));

    if ( nargin == 1 )
        set(gcf,'position', [150 341 830 559]);
    else
        set(gcf,'position', [150 341 1000 559]);
    end

    inp = input(sprintf('#%.0f of %.0f: Keep or reject? <enter> = reject, 1 = keep: ',...
        i,length(ccpairs)));
    if ( ~isempty(inp) )
        index = [index i];
        save('temp-ccpairs-select', 'index');
    end


end

return;




