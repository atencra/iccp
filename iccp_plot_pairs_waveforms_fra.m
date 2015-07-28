function vs_plot_pairs_waveforms_fra(fra)
%plot_strf - Plots STRF data obtained using the Michigan 16
%   channel silicon probe.
%
%   plot_strf(strf, trigger, flim, tlim)
%
%   strf is a struct array holding the receptive field data
%
%   trigger is a vector of trigger times associated with the
%     ripple stimlui
%
%   flim and tlim are optional input arguments.
%
%   flim is an optional 1x2 vector, holding the frequency range 
%      over which the strfs are plotted. If flim is not input then 
%      the strfs are plotted over the full frequency range.
%
%   tlim is the same as flim, except it is over the time axis.
%
%   caa 6/28/02

if ( nargin == 4 )
   flim = [];
   tlim = [];
end

if ( nargin == 5 )
   tlim = [];
end

close all;
%3,15
for i = 9:9 %length(ccpairs)
i
   figure;

   exp = ccpairs(i).exp;
   site = ccpairs(i).site;
   chan = ccpairs(i).chan;
   model1 = ccpairs(i).model1;
   model2 = ccpairs(i).model2;
   position = ccpairs(i).position;
   stim = ccpairs(i).stim;

   pd = ccpairs(i).peakdelay;
   hw = ccpairs(i).halfwidth;
   ccc = ccpairs(i).ccc;
   significant = ccpairs(i).significant;

   conf_limit = ccpairs(i).conf_limit;
   delay = ccpairs(i).delay;
   q12 = ccpairs(i).q12;

   nsp1 = length(ccpairs(i).spiketimes1);
   nsp2 = length(ccpairs(i).spiketimes2);


   % Get STA/nonlinearities that match the neurons in the current
   % ccpairs(i)
   [fio1, fio2] = get_sta_fio_from_ccpairs_neuron(ccpairs(i), fioFit);
   [strf1, strf2] = get_strf_from_ccpairs_neuron(ccpairs(i), strf);


   % Get spike waveforms
   % ----------------------------------------
   wavetime1 = ccpairs(i).meanwaveform1(:,1);
   waveform1 = ccpairs(i).meanwaveform1(:,2);
   waveform1= waveform1 ./ 2048;
   waveform1 = smooth3p(waveform1);

   wavetime2 = ccpairs(i).meanwaveform2(:,1);
   waveform2 = ccpairs(i).meanwaveform2(:,2);
   waveform2= waveform2 ./ 2048;
   waveform2 = smooth3p(waveform2);

   % Normalize the waveforms to [-1, 1]
   max1 = max([max(abs(waveform1)) max(abs(min(waveform1)))]);
   max2 = max([max(abs(waveform2)) max(abs(min(waveform2)))]);
   maxmax = max([max1 max2]);
   waveform1 = waveform1 ./ maxmax;
   waveform2 = waveform2 ./ maxmax;

   % Get STAs
   sta1 = zeros(size(fio1.sta{1}));
   for j = 1:length(fio1.sta)
      sta1 = sta1 + fio1.sta{j};
   end % (for j)

   sta2 = zeros(size(fio2.sta{1}));
   for j = 1:length(fio2.sta)
      sta2 = sta2 + fio2.sta{j};
   end % (for j)

   % Get nonlinearities
   dt = fio1.dt;
   x_sta1 = fio1.x_sta;
   fx_sta1 = fio1.fx_sta ./ dt;
   xFit_sta1 = fio1.xFit_sta;
   fxFit_sta1 = fio1.fxFit_sta ./ dt;

   x_sta2 = fio2.x_sta;
   fx_sta2 = fio2.fx_sta ./ dt;
   xFit_sta2 = fio2.xFit_sta;
   fxFit_sta2 = fio2.fxFit_sta ./ dt;

   nmse1 = fio1.nmse_sta;
   nmse2 = fio2.nmse_sta;

   theta1 = fio1.fitParams_sta(2);
   sigma1 = fio1.fitParams_sta(3);

   theta2 = fio2.fitParams_sta(2);
   sigma2 = fio2.fitParams_sta(3);


 
   % Plot the waveforms
   % -------------------------------
   hs1 = subplot(3,2,1);
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
   title(sprintf('#%.0f of %.0f: %s s%.0fc%.0fm%.0f-%.0f', ...
      i, length(ccpairs), exp, site, chan, model1, model2));



   % Plot the cross-covariance
   % -------------------------------
   hs2 = subplot(3,2,2);
   upper95qab = conf_limit;
   lower95qab = -conf_limit;
   hold on;
   ymin = min([min(q12) lower95qab]);
   ymax = max([max(q12) upper95qab]);
   range = ymax-ymin;
   bar(delay, q12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
   plot([pd pd], ylim, 'r-');
   plot([pd-hw/2 pd-hw/2], ylim, 'b-');
   plot([pd+hw/2 pd+hw/2], ylim, 'b-');
   plot([0 0], [ymin-0.1*range ymax+0.1*range], 'k-');
   plot([min(delay) max(delay)], [0 0], 'k-');
   plot([min(delay) max(delay)], [upper95qab upper95qab], 'r-');
   plot([min(delay) max(delay)], [lower95qab lower95qab], 'r-');
   xlim([min(delay) max(delay)]);
   xlim([-25 25]);
   xtick = -25:5:25;
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   xlabel('Delay (ms)');
   ylabel('Cross Cov');
   ylim([ymin-0.1*range ymax+0.1*range]);
   tickpref;
   title(sprintf('sig = %.0f, pd = %.1f, hw = %.1f', significant, pd, hw));



   subplot(3,2,3);
   plot_strf_single(strf1,trigger,flim,tlim);

   subplot(3,2,4);
   plot_strf_single(strf2,trigger,flim,tlim);



   subplot(3,2,5);
   hold on;
   plot(xFit_sta1, fxFit_sta1, '-', 'color', 0.6*ones(1,3), 'linewidth', 1.5);
   plot(x_sta1, fx_sta1, 'ko', 'markerfacecolor', 'k', 'markersize', 2);
   xlim([min(x_sta1) max(x_sta1)]);
   ylim([0 max(fx_sta1)]);
   box off;
   tickpref;
%    set(gca,'xtick', [-1:2], 'xticklabel', [-1:2]);
   xlabel('Projection (SD)');
%    set(gca,'yticklabel', [-1:1], 'ytick', [-1:1]);
   ylabel('Firing Rate (sp/s)');
   title(sprintf('%s; nmse=%.2f,t=%.2f,s=%.2f', stim, nmse1, theta1, sigma1));


   subplot(3,2,6);
   hold on;
   plot(xFit_sta2, fxFit_sta2, '-', 'color', 0.6*ones(1,3), 'linewidth', 1.5);
   plot(x_sta2, fx_sta2, 'ko', 'markerfacecolor', 'k', 'markersize', 2);
   xlim([min(x_sta2) max(x_sta2)]);
   ylim([0 max(fx_sta2)]);
   box off;
   tickpref;
%    set(gca,'xtick', [-1:2], 'xticklabel', [-1:2]);
   xlabel('Projection (SD)');
%    set(gca,'yticklabel', [-1:1], 'ytick', [-1:1]);
   ylabel('Firing Rate (sp/s)');
   title(sprintf('nmse=%.2f,t=%.2f,s=%.2f', nmse2, theta2, sigma2));

   set(gcf,'position',[450 30 438 571]);

end

return;





function [fio1, fio2] = get_sta_fio_from_ccpairs_neuron(ccpairs, fioFit)

exp = ccpairs.exp;
site = ccpairs.site;
chan = ccpairs.chan;
model1 = ccpairs.model1;
model2 = ccpairs.model2;
position = ccpairs.position;
stim = ccpairs.stim;

found = 0;
i = 1;
while ( i <= length(fioFit) && ~found )
   exp_fio = fioFit(i).exp;
   site_fio = fioFit(i).site;
   chan_fio = fioFit(i).chan;
   model_fio = fioFit(i).model;
   stim_fio = fioFit(i).stim;

% [strcmp(exp, exp_fio)]
% [site site_fio]
% [chan chan_fio]
% [sum(model1 - model_fio)]
% [strcmp(stim, stim_fio)]

   if ( strcmp(exp, exp_fio) && site == site_fio && chan == chan_fio && ...
      sum(model1 - model_fio) == 0 && strcmp(stim, stim_fio) )
      fio1 = fioFit(i);
      found = 1;
   end

   i = i + 1;
end % (while)


found = 0;
i = 1;
while ( i <= length(fioFit) && ~found )
   exp_fio = fioFit(i).exp;
   site_fio = fioFit(i).site;
   chan_fio = fioFit(i).chan;
   model_fio = fioFit(i).model;
   stim_fio = fioFit(i).stim;

   if ( strcmp(exp, exp_fio) && site == site_fio && chan == chan_fio && ...
      sum(model2 - model_fio) == 0 && strcmp(stim, stim_fio) )
      fio2 = fioFit(i);
      found = 1;
   end
   i = i + 1;
end % (while)

return;




function [strf1, strf2] = get_strf_from_ccpairs_neuron(ccpairs, strf)

exp = ccpairs.exp;
site = ccpairs.site;
chan = ccpairs.chan;
model1 = ccpairs.model1;
model2 = ccpairs.model2;
position = ccpairs.position;
stim = ccpairs.stim;

found = 0;
i = 1;
while ( i <= length(strf) && ~found )
   exp_strf = strf(i).exp;
   site_strf = strf(i).site;
   chan_strf = strf(i).chan;
   model_strf = strf(i).model;
   stim_strf = strf(i).stim;

   if ( strcmp(exp, exp_strf) && site == site_strf && chan == chan_strf && ...
      sum(model1 - model_strf) == 0 && strcmp(stim, stim_strf) )
      strf1 = strf(i);
      found = 1;
   end

   i = i + 1;
end % (while)


found = 0;
i = 1;
while ( i <= length(strf) && ~found )
   exp_strf = strf(i).exp;
   site_strf = strf(i).site;
   chan_strf = strf(i).chan;
   model_strf = strf(i).model;
   stim_strf = strf(i).stim;

   if ( strcmp(exp, exp_strf) && site == site_strf && chan == chan_strf && ...
      sum(model2 - model_strf) == 0 && strcmp(stim, stim_strf) )
      strf2 = strf(i);
      found = 1;
   end
   i = i + 1;
end % (while)

return;



