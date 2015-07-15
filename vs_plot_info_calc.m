function iresp = vs_plot_info_calc(resp)
% NEEDS TO BE MODULARIZED
%
% info_calculation: Calculates information values and nmse values from an
% STRF response file. (saved as resp_final)
%
% [info_nmse, sigma_nmse, nmse] = info_calculation(resp, graphs, popgraphs)
%
% Inputs:
% -----------------------------
% 1) resp: STRF response file, saved as -*final-resp. 
% 2) graphs: boolean input. 1 if graphs are desired, 0 if not.
% 3) popgraphs: boolean input. 1 if population graphs are desired, 0 if
%       not.
%
% Outputs:
%   ----------------------------
%   IF GRAPHS is 1 (plots are for each individual neuron):
%
%   1) Generates a plot of Information vs. all the dt values.
%   2) Generates a PSTH, raster plot, and image of the stimulus 
%             spectrogram, using the value of dt where the I vs dt graph 
%             becomes nonlinear.
%   3) Plots all of its information values, divided into intervals of 
%             width [[0.1 0.5 1 2 3 4 6 12], against log(1/those T values).
%   4) Plots sigma against 1/sqrt(those T values).
%
%   IF POPGRAPHS is 1 (plots are for all the neurons)
%     1) Plots the population histogram of normalized mean square error, 
%             information, and picked dt values.
%
%  Finally, performs pairwise analysis on information, sigma values, and
%  generates scatterplots of those.

iresp = resp;

total_dur = 12;
dt = linspace(0.0001,0.050,100);

for k = 1:length(resp)

   trialspikes = resp(k).trialspikes; % cell array of spike responses

   % Get preliminary psth to see if there are any spikes
   [p, rt] = vs_psth_firing_rate(trialspikes, dt(end), total_dur);

   if ( length(find(rt>0)) > 0 ) % There were spikes

      [idt, dt] = vs_info_time_bin( trialspikes, dt, total_dur ); % information for time bin sizes

      iresp(k).dt = dt;
      iresp(k).idt = idt;

      [dtoptim] = vs_info_pick_dt(idt, dt); % optimal time bin size, in seconds
      iresp(k).dtoptim = dtoptim;

      [p, rt] = vs_psth_firing_rate(trialspikes, dtoptim, total_dur); % psth, firing rate

      rbar = mean(rt); % average firing rate

      frac_dur = total_dur * [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.85 0.9 0.925 0.95 0.975 1];
      iresp(k).total_dur = total_dur;
      iresp(k).frac_dur = frac_dur;

      ifracdur = cell(1,length(frac_dur));

      for j = 1:length(frac_dur)

         sampdur = frac_dur(j); % reduced duration, in seconds
         nvals = floor(sampdur ./ dtoptim); % number of bins in reduced duration
         maxval = length(rt)-1; % maximum number of time bins

         for ii = 1:100
            index = randperm(maxval, nvals); % choose samples from maximum number of bins
            rt_subset = rt(index);
            index = find( rt_subset > 0 );
            irt(ii) = dtoptim * 1./sampdur .* ...
               sum( (rt_subset(index)./rbar).*log2( rt_subset(index)./rbar) );
         end % (for ii)

         ifracdur{j} = irt;

      end % (for j)

      iresp(k).ifracdur = ifracdur; 

      figure;
      subplot(3,1,1);

      ifracdurmn = zeros(1,length(ifracdur));
      ifracdurst = zeros(1,length(ifracdur));
      for ii = 1:length(ifracdur)
         ivals = ifracdur{ii};
         ifracdurmn(ii) = mean(ivals);
         ifracdursd(ii) = std(ivals);
         hold on;
         plot(ones(size(ivals)).*1./frac_dur(ii), ivals, 'ko');
      end

      iresp(k).ifracdurmn = ifracdurmn;
      iresp(k).ifracdursd = ifracdursd;

      set(gca,'xscale', 'log');
      xlim([min(1./frac_dur) max(1./frac_dur)]);
      tickpref;
      xlabel('1/T');
      ylabel('Info');

      subplot(3,1,2);
      plot(1./frac_dur, ifracdurmn, 'ko');
      set(gca,'xscale', 'log');
      xlim([min(1./frac_dur) max(1./frac_dur)]);
      ylim([0 1.1*max(ifracdurmn)]);
      tickpref;
      xlabel('1/T');
      ylabel('Info');

      subplot(3,1,3);
      x = 1./sqrt(frac_dur);
      y = ifracdursd;
      beta = polyfit(x, y, 1);
      xfit = linspace(min(x), max(x), 500);
      yfit = polyval(beta, xfit);
      hold on;
      plot(x, y, 'ko');
      plot(xfit, yfit, 'k-');
      xlim([sqrt(min(1./frac_dur)) sqrt(max(1./frac_dur))]);
      ylim([0 1.1*max(ifracdursd)]);
      tickpref;
      xlabel('1/sqrt(T)');
      ylabel('Info SD');

      set(gcf,'position', [245   193   415   799]);

   end % (if)

end % (for k)

return;






