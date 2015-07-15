function [r] = vs_plot_psth(resp, dt, T)

% plot_psth Generates raster plot, PSTH, and stimulus spectrogram for a
% given response file and dt value

% [r] = plot_psth(resp,dt,graphs)

% Inputs:
% ------------------------------
% 1) resp: STRF response data. Saved as resp_final, or input from
%       info_calculation. Length is 1 if called from info_calculation.
% 2) dt: dt value used as a bin width
% 3) graphs: boolean input. Outputs the PSTH, raster, and stimulus
%         spectrogram plots if 1, does not do so if 0.

% Outputs:
% ------------------------------
% 1) r: Cell array containing the number of spikes within each bin of width
%       dt, used to calculate information values.
% 2) Plot of the stimulus spectrogram.
% 3) Raster plot of all the trialspikes values in the resp
%       file.
% 4) Population histogram of all spiketimes.


L = length(resp);
edges = 0:dt:T;
r = cell(1,L);

if graphs
    sprfile = 'C:\Users\Victor\Documents\Research\20040217\stimuli_20040217\dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-12sec.spr';
    [stimulus] = read_entire_spr_file(sprfile);
    len = length(stimulus);
    timeaxis = 0:4009:23280;
end


for i = 1:L
    all_spiketimes = [];
    
    if graphs
        figure;
        subplot(3,4,[1 4])
        imagesc(stimulus)
        set(gca,'XTick',timeaxis,'XTickLabel',timeaxis)
        suptitle(sprintf('PSTH and Spectrogram, 4009 samples = 2 seconds, dt = %.3f seconds',dt))
        box off
        tickpref
    end
    
    %     Current File Names, NEED TO CONVERT
    for j = 1:length(resp(i).trialspikes)
        spikes = resp(i).trialspikes{j}./1000;
        all_spiketimes = [all_spiketimes spikes];
        
        if graphs
            subplot(3,4,[5 8])
            hold on
            plot(spikes,j*ones(length(spikes)),'k.','MarkerSize',4)
        end
         box off
        tickpref
    end
    
    % Standard File Names
    %  for j = 1:length(data(i).resp)
    %         spikes = data(i).resp{j}./1000;
    %         all_spiketimes = [all_spiketimes spikes];
    %
    %         subplot(4,2,[1 4])
    %         hold on
    %         plot(spikes,j*ones(length(spikes)),'k.')
    %
    %  end
    
    n = histc(all_spiketimes, edges);
    
    if graphs
        subplot(3,4,[9 12])
        hb = bar(edges,n,'histc');
        set(hb, 'FaceColor', [0.6 0.6 0.6]);
        xlim([0 12])
        box off
        tickpref
        pause;
        close
    end
    
    r{i} = n; % number of spikes within each bin of width dt
     
end