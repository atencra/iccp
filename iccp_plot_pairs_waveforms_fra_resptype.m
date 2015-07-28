function vs_plot_pairs_waveforms_fra_resptype(spk, fra, params, rt)
% vs_plot_pairs_fra_pairs_resptype :
% --------------------------------------------------------
% Reads through files in a directory to find struct arrays
% holding the pairs data. The files have names in the following
% form:
%        *-fra-params.mat
%
% params : struct array holding fra parameters, such
% as characteristic frequency, bandwidth, latency and firing rate.
%
% The function finds the neurons that were recorded from the 
% same channel, and then compares the parameters
% for these neurons.
%
% The analysis tells us the variability of the parameters in
% a locally confined area of the cortex.
%
% caa 6/29/12
% 

chan = [spk.chan];
[temp, index] = sort(chan);
spk = spk(index);

chan = [fra.chan];
[temp, index] = sort(chan);
fra = fra(index);

chan = [params.chan];
[temp, index] = sort(chan);
params = params(index);

chan = [rt.chan];
[temp, index] = sort(chan);
rt = rt(index);


chan = [fra.chan];
chan_unique = unique(chan);

cmb_total = [];
for i = 1:length(chan_unique)
   index = find( chan_unique(i) == chan );
   cmb = nchoosek(index, 2); % determine all possible pairwise combinations
   [nr, nc] = size(cmb);
   cmb_total = [cmb_total; cmb];
end

close all;
%2003-11-24 site9: 17, 19, 20, 26, 38, 31, 37, 
%2003-11-24 site10: 32, 34,
%2003-11-24 site11: 3, 7, 9, 13, 22, 
%2003-11-24 site14: 1, 5, 20, 42,
for i = 42:42 %1:size(cmb_total,1)

   figure;

   index1 = cmb_total(i,1);
   index2 = cmb_total(i,2);

   % Get data for the first unit in the pair
   % ----------------------------------------
   fra1 = fra(index1);
   [f1, a1, tc1] = freq_resp_area_matrix(fra1, [0 55]);
   spk1 = spk(index1);
   params1 = params(index1);
   rt1 = rt(index1);
   exp1 = fra1.exp; % experiment
   site1 = fra1.site; % penetration number
   stim1 = fra1.stim; % type of stimulus
   depth1 = fra1.depth; % depth of probe tip
   chan1 = fra1.chan;
   model1 = fra1.model;
   position1 = fra1.position;
   nreps1 = fra1.nreps;
   psi1 = rt1.ne_nb_ne_total;



   % Get data for the second unit in the pair
   % ----------------------------------------
   fra2 = fra(index2);
   [f2, a2, tc2] = freq_resp_area_matrix(fra2, [0 55]);
   spk2 = spk(index2);
   params2 = params(index2);
   rt2 = rt(index2);
   exp2 = fra2.exp; % experiment
   site2 = fra2.site; % penetration number
   stim2 = fra2.stim; % type of stimulus
   depth2 = fra2.depth; % depth of probe tip
   chan2 = fra2.chan;
   model2 = fra2.model;
   position2 = fra2.position;
   nreps2 = fra2.nreps;
   psi2 = rt2.ne_nb_ne_total;


   % Get some preliminary data.
   amp = fra1.amp; % max sound amplitude
   freq = fra1.freq; % vector of frequencies
   aunq = unique(amp);
   aunq = aunq(:)';
   funq = unique(freq);
   funq = funq(:)';


   % Data for nice tick labels for the plots.
   ttick = [1 5 10 15];
   temp = unique(amp);
   tlabel = temp(ttick);
   ftick = [1 11 22 33 45];
   ftick = [1 23 45];
   temp = unique(freq);
   flabel = round(10*temp(ftick))/10;


   % Get spike waveforms
   % ----------------------------------------
   wavetime1 = spk(index1).meanwaveform(:,1);
   waveform1 = spk(index1).meanwaveform(:,2);
   waveform1= waveform1 ./ 2048;
   waveform1 = smooth3p(waveform1);

   wavetime2 = spk(index2).meanwaveform(:,1);
   waveform2 = spk(index2).meanwaveform(:,2);
   waveform2= waveform2 ./ 2048;
   waveform2 = smooth3p(waveform2);


   % Normalize the waveforms to [-1, 1]
   max1 = max([max(abs(waveform1)) max(abs(min(waveform1)))]);
   max2 = max([max(abs(waveform2)) max(abs(min(waveform2)))]);
   maxmax = max([max1 max2]);
   waveform1 = waveform1 ./ maxmax;
   waveform2 = waveform2 ./ maxmax;



   % Plot the waveforms
   % -------------------------------
   subplot(2, 3, 1);
   hold on;
   plot(wavetime1, waveform1, 'k-', 'linewidth', 1.5);
   plot(wavetime2, waveform2, '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5);
   legend('1', '2');
   xlim([-1 2]);
   ylim([-1.05 1.05]);
   box off;
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca,'xtick', [-1:2], 'xticklabel', [-1:2]);
   xlabel('Time (ms)');
   set(gca,'yticklabel', [-1:1], 'ytick', [-1:1]);
   ylabel('Norm. Amp.');
   title(sprintf('%s', exp1));


   % Plot the fra for unit 1
   % -------------------------------
   subplot(2, 3, 2);
   imagesc(tc1);
   axis xy;
   tickpref;
   colormap(flipud(colormap(gray)));
   set(gca,'xtick', ftick, 'xticklabel', flabel);
   xlabel('Frequency (kHz)');
   set(gca,'ytick', ttick, 'yticklabel', tlabel);
   ylabel('dB SPL');
   for j = 1:length(model1)
      m1(j) = num2str(model1(j));
   end % (for j)
   title(sprintf('s%.0f c%.0f m%s', site1, chan1, m1));


   % Plot the fra for unit 2
   % -------------------------------
   subplot(2, 3, 5);
   imagesc(tc2);
   axis xy;
   tickpref;
   colormap(flipud(colormap(gray)));
   set(gca,'xtick', ftick, 'xticklabel', flabel);
   xlabel('Frequency (kHz)');
   set(gca,'ytick', ttick, 'yticklabel', tlabel);
   ylabel('dB SPL');
   for j = 1:length(model2)
      m2(j) = num2str(model2(j));
   end % (for j)
   title(sprintf('s%.0f c%.0f m%s', site2, chan2, m2));

   cmap = brewmaps('reds', 25);
   cmap = flipud(cmap);
   cmap = [1 1 1; cmap];
   colormap(cmap);




   psth1 = [];
   for ii = 1:length(aunq)
      for iii = 1:length(funq)
         i1 = find( freq==funq(iii) & amp==aunq(ii) );
         for iv = 1:length(i1)
            raster{(iii-1)*nreps1 + iv} = fra1.resp{i1(iv)};
            psth1 = [psth1 fra1.resp{i1(iv)}];
         end % (for j)
      end % (for i)
   end % (for m)

   subplot(2,3,3);
   hold on;
   edges = 0:5:300;
   count1 = histc(psth1,edges);
   fr1  = count1 / 0.005 / nreps1 / length(aunq) / length(funq);
   p = patch([1 1 50 50],[0 max(fr1) max(fr1) 0],[0.9 0.9 0.9]);
   set(p,'linestyle', 'none');
   hb = bar(edges, fr1, 'histc');
   set(hb, 'facecolor', [0.3 0.3 0.3]);
   set(gca,'xtick', 0:150:300, 'xticklabel', 0:150:300);
   tickpref;
   xlim([0 300]);
   ylim([0 1.05*max(fr1)]);
   xlabel('Time (ms)');
   ylabel('Firing Rate (sp/s)');
   title(sprintf('PSI = %.3f',psi1));


   psth2 = [];
   for ii = 1:length(aunq)
      for iii = 1:length(funq)
         i1 = find( freq==funq(iii) & amp==aunq(ii) );
         for iv = 1:length(i1)
            raster{(iii-1)*nreps2 + iv} = fra2.resp{i1(iv)};
            psth2 = [psth2 fra2.resp{i1(iv)}];
         end % (for j)
      end % (for i)
   end % (for m)

   subplot(2,3,6);
   hold on;
   edges = 0:5:300;
   count2 = histc(psth2,edges);
   fr2  = count2 / 0.005 / nreps2 / length(aunq) / length(funq);
   p = patch([1 1 50 50],[0 max(fr2) max(fr2) 0],[0.9 0.9 0.9]);
   set(p,'linestyle', 'none');
   hb = bar(edges, fr2, 'histc');
   set(hb, 'facecolor', [0.3 0.3 0.3]);
   set(gca,'xtick', 0:150:300, 'xticklabel', 0:150:300);
   tickpref;
   xlim([0 300]);
   ylim([0 1.05*max(fr2)]);
   xlabel('Time (ms)');
   ylabel('Firing Rate (sp/s)');
   title(sprintf('PSI = %.3f',psi2));


   set(gcf,'position', [278 626 438 246]);
   print_mfilename(mfilename);

end

return;




% % get the latency data
% % -------------------------------------------------------------
% lat = [];
% for i = 1:length(fra(n).resp)
% temp = fra(n).resp{i};
% index = find(temp>0 & temp<50);
% lat = [lat temp(index)];
% end % (for i)
% fig_latency = figure;
% hist(lat,50); % 50 bins gives about 1 msec resolution
% set(gcf,'position', [182   525   560   420]);
% [latency, dummay] = ginput(1);
% close(fig_latency);
% latency



function [position, nb_nb_ne, ne_nb_ne] = get_pairs_fra_resptype

position = [];
nb_nb_ne = [];
ne_nb_ne = [];

d = dir('*-fracmb-pairs-resptype.mat');

for n = 1:length(d)

    filename = d(n).name;
    s = load(filename, 'rt');
    rt = s.rt;
    
    chan = [rt.chan];
    chan_unique = unique(chan);
    
    for i = 1:length(chan_unique)
        
        index = find( chan_unique(i) == chan );
        
        if length(index)>1
            
            cmb = nchoosek(index, 2); % determine all possible pairwise combinations
            
            [nr, nc] = size(cmb);
            
            for j = 1:nr
                
                index1 = cmb(j,1);
                index2 = cmb(j,2);
                
                position = [position; rt(index1).position];
                
                nb_nb_ne = [nb_nb_ne; rt(index1).nb_nb_ne_total rt(index2).nb_nb_ne_total];
                ne_nb_ne = [ne_nb_ne; rt(index1).ne_nb_ne_total rt(index2).ne_nb_ne_total];
                
            end % (for j)
            
        end % (if/else)
        
    end % (for i)

end % (for n)
return;



function vs_plot_pairs_fra_resptype_scatter_hist(nb_nb_ne, ne_nb_ne)

rtdata{1} = nb_nb_ne;
rtdata{2} = ne_nb_ne;

nb_nb_ne_diff = abs( nb_nb_ne(:,1) - nb_nb_ne(:,2) ); 
ne_nb_ne_diff = abs( ne_nb_ne(:,1) - ne_nb_ne(:,2) ); 

mn_nb_nb_ne_diff = mean( nb_nb_ne_diff );
md_nb_nb_ne_diff = median( nb_nb_ne_diff );
sd_nb_nb_ne_diff = std( nb_nb_ne_diff );

mn_ne_nb_ne_diff = mean( ne_nb_ne_diff );
md_ne_nb_ne_diff = median( ne_nb_ne_diff );
sd_ne_nb_ne_diff = std( ne_nb_ne_diff );

rtdiff{1} = nb_nb_ne_diff;
rtdiff{2} = ne_nb_ne_diff;

rtdiff_mn{1} = mn_nb_nb_ne_diff;
rtdiff_mn{2} = mn_ne_nb_ne_diff;

rtdiff_md{1} = md_nb_nb_ne_diff;
rtdiff_md{2} = md_ne_nb_ne_diff;

rtdiff_sd{1} = sd_nb_nb_ne_diff;
rtdiff_sd{2} = sd_ne_nb_ne_diff;




minmin = 0;
maxmax = 1;

figure;
edges = linspace(0,1,21);

label = {'NB-NB-NE', 'NE-NB-NE'};
tick = [0:0.25:1];
nrows = length(label);

for i = 1:nrows

   subplot(nrows,2,(i-1)*2+1);
   hold on
   xlim([minmin maxmax]);
   ylim([minmin maxmax]);
   plot(xlim, ylim, 'k-');
   [larger, smaller] = vs_largersmaller(rtdata{i}(:,1),rtdata{i}(:,2));
   scatter(larger,smaller,20,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.6 0.6 0.6])
   xlim([minmin maxmax])
   ylim([minmin maxmax])
   tickpref
   box off
   set(gca,'xtick', tick, 'xticklabel', tick);
   set(gca,'ytick', tick, 'yticklabel', tick);
   ylabel(sprintf('%s Neuron 2',label{i}))
   xlabel(sprintf('%s Neuron 1',label{i}))

   subplot(nrows,2,(i)*2);
   n = histc(rtdiff{i}, edges);
   n = n./sum(n);
   hb = bar(edges,n,'histc');
   set(hb, 'FaceColor', [0.6 0.6 0.6])
   xlabel('Response Type Difference')
   xlim([min(edges) max(edges)]);
   ylim([0 0.4]);
   set(gca,'xtick', tick, 'xticklabel', tick);
   ytick = 0:0.1:0.7;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   tickpref
   box off



   [r,p] = corrcoef(larger, smaller );
   fprintf('\n');
   fprintf('%s\n',label{i});
   fprintf('r = %.4f, p = %.4f\n',r(2), p(2) );

   fprintf('\n');
   fprintf('%s Difference\n', label{i});
   fprintf('MN = %.4f, MD = %.4f, SD = %.4f, SE = %.4f, n=%.0f\n',...
      rtdiff_mn{i}, rtdiff_md{i},rtdiff_sd{i},rtdiff_sd{i}./sqrt(length(rtdiff_sd{i})), length(rtdata{i}(:,1)) );

end % (for i)

set(gcf,'Position',[100 100 438 700])
print_mfilename(mfilename);

return;


















function plot_strf_spk_pairs(strf, trigger, spk, flim, tlim)
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
%   spk : struct array holding spike waveforms and spike times
%   Must be the same length as strf, where each element of spk
%   corresponds to an element in strf.
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
%
% plot_strf_spk_pairs(strf, trigger, spk, flim, tlim)

if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   flim = [];
   tlim = [];
end

if ( nargin == 4 )
   tlim = [];
end

time = 1000 .* strf(1).taxis; % change to ms
freq = strf(1).faxis;
depth = strf(1).depth;
site = strf(1).site;
exp = strf(1).exp;
stim = strf(1).stim;


% Now get some time tick marks
if ( isempty(tlim) )
   t = round(min(time)):abs(round(min(time))):max(time);
   ttick = [];
   for i = round(min(time)):abs(round(min(time))):max(time)
      temp = abs(i+0.01-time);
      ttick = [ttick find(temp == min(temp))];
   end
else

   dttick = ( abs(tlim(1)) + abs(tlim(2)) ) / 5;
   t = round(tlim(1)):dttick:round(tlim(2));

   index_tlim_min = find( time >= tlim(1) );
   index_tlim_min = min(index_tlim_min);

   index_tlim_max = find( time <= tlim(2) );
   index_tlim_max = max(index_tlim_max);

   ttick = [];
   for i = round(tlim(1)):dttick:round(tlim(2))
      temp = abs(i+0.01-time);
      ttick = [ttick find(temp == min(temp))];
   end
   ttick = ttick - (index_tlim_min - 1);
end


% Now we get some frequency tick marks
if ( isempty(flim) )
   nocts = floor(log2(max(freq)/min(freq)));
   f = round(min(freq).*2.^(0:2:nocts));
   ftick = [];
   for i = 0:2:nocts
      temp = abs(i+0.001-log2(freq/min(freq)));
      ftick = [ftick find(temp == min(temp))];
   end
else

   nocts = floor(log2(max(flim(2))/min(flim(1))));
   f = round(min(flim(1)).*2.^(0:1:nocts));
   ftick = [];
   for i = 0:1:nocts
      temp = abs(i+0.001-log2(freq/min(flim(1))));
      ftick = [ftick find(temp == min(temp))];
   end
   index_flim_min = find( freq >= flim(1) );
   index_flim_max = find( freq <= flim(2) );
   index_flim_min = min(index_flim_min);
   index_flim_max = max(index_flim_max);
   ftick = ftick - (index_flim_min - 1);
end


fs = strf(1).fs; % sampling rate of A/D system
mdb = strf(1).mdb;
sigma2 = (mdb)^2 / 8; % variance of dynamic moving ripple
sigma = sqrt(sigma2); % std of dynamic moving ripple
dur = ( trigger(end)-trigger(1) ) / fs; % total duration of ripple, in sec

ntot = 1;
len = length(strf);
fignum = 1;
figure;

gaussian = fspecial('gaussian', [20 20], 7);
% imagesc(gaussian)
% pause

p = 0.002;

close all;

chan = [strf.chan];
chan_unique = unique(chan);

cmb_total = [];
for i = 1:length(chan_unique)
   index = find( chan_unique(i) == chan )
   cmb = nchoosek(index, 2); % determine all possible pairwise combinations
   [nr, nc] = size(cmb);
   cmb_total = [cmb_total; cmb];
end

figure_row = 1;

for i = 1:size(cmb_total,1)

      index1 = cmb_total(i,1);
      index2 = cmb_total(i,2);

      % Get data for the first unit in the pair
      % ----------------------------------------
      chan1 = strf(index1).chan;
      model1 = strf(index1).model;
      position1 = strf(index1).position;

      rf1 = double(strf(index1).rfcontra);
      n01 = strf(index1).n0contra;
      w01 = strf(index1).w0contra;
      [rfsig1] = significant_strf(rf1, p, n01, mdb, dur);

      % Get data for the second unit in the pair
      % ----------------------------------------
      chan2 = strf(index2).chan;
      model2 = strf(index2).model;
      position2 = strf(index2).position;

      rf2 = double(strf(index2).rfcontra);
      n02 = strf(index2).n0contra;
      w02 = strf(index2).w0contra;
      [rfsig2] = significant_strf(rf2, p, n02, mdb, dur);


      % Shift strfs for proper time axis position
      % ----------------------------------------
%      rfsig1 = circshift(rfsig1, [0 120]);
%      rfsig2 = circshift(rfsig2, [0 120]);



      % Cut axes to fit in desired range
      % ----------------------------------------
      if ( ~isempty(flim) )
         rfsig1 = rfsig1(index_flim_min:index_flim_max, :);
      end

      if ( ~isempty(tlim) )
         rfsig1 = rfsig1(:, index_tlim_min:index_tlim_max);
      end

      if isempty(rfsig1)
         rfsig1 = zeros(length(freq), length(time));
      end

      if ( ~isempty(flim) )
         rfsig2 = rfsig2(index_flim_min:index_flim_max, :);
      end

      if ( ~isempty(tlim) )
         rfsig2 = rfsig2(:, index_tlim_min:index_tlim_max);
      end

      if isempty(rfsig2)
         rfsig2 = zeros(length(freq), length(time));
      end




      % Get spike waveforms
      % ----------------------------------------
      wavetime1 = spk(index1).meanwaveform(:,1);
      waveform1 = spk(index1).meanwaveform(:,2);
      waveform1= waveform1 ./ 2048;

      wavetime2 = spk(index2).meanwaveform(:,1);
      waveform2 = spk(index2).meanwaveform(:,2);
      waveform2= waveform2 ./ 2048;


      % Normalize the waveforms to [-1, 1]
      max1 = max([max(abs(waveform1)) max(abs(min(waveform1)))]);
      max2 = max([max(abs(waveform2)) max(abs(min(waveform2)))]);
      maxmax = max([max1 max2]);
      waveform1 = waveform1 ./ maxmax;
      waveform2 = waveform2 ./ maxmax;

 

      % Plot the waveforms
      % -------------------------------
      hs1 = subplot(5,3, 3*(figure_row-1) + 1);
      hold on;
      plot(wavetime1, waveform1, 'k-', 'linewidth', 1.5);
      plot(wavetime2, waveform2, '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5);
      if ( figure_row == 1 )
         legend('1', '2');
      end
      xlim([-1 2]);
      ylim([-1.05 1.05]);
      box off;
      set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);

      set(gca,'xtick', [-1:2], 'xticklabel', [-1:2]);
      if ( figure_row == 5 || i == size(cmb_total,1) )
         xlabel('Time (ms)');
      end
      set(gca,'yticklabel', [-1:1], 'ytick', [-1:1]);
      ylabel('Norm. Amp.');
   
   
      % Plot the strf for unit 1
      % -------------------------------
      subplot(5,3, 3*(figure_row-1) + 2);
      imagesc(rfsig1);
      box off;
      set(gca,'ydir', 'normal');
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      minmin = min(min(rfsig1));
      maxmax = max(max(rfsig1));
      boundary1 = max([abs(minmin) abs(maxmax)]);
      set(gca, 'clim', [-1.05*boundary1-eps 1.05*boundary1+eps]);
      set(gca,'xtick', ttick, 'xticklabel',t);
      if ( figure_row == 5 || i == size(cmb_total,1) )
         xlabel('Time (ms)');
      end
      set(gca,'ytick', ftick, 'yticklabel', f./1000);
      title(sprintf('1: %.0f-%.0f %.0fum %.2f Hz', chan1, model1, position1, w01));
      ylabel(sprintf('Frequency (kHz)'));


      % Plot the strf for unit 2
      % -------------------------------
      subplot(5,3, 3*(figure_row-1) + 3);
      imagesc(rfsig2);
      box off;
      set(gca,'ydir', 'normal');
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      minmin = min(min(rfsig2));
      maxmax = max(max(rfsig2));
      boundary2 = max([abs(minmin) abs(maxmax)]);
      set(gca, 'clim', [-1.05*boundary2-eps 1.05*boundary2+eps]);
      set(gca,'xtick', ttick, 'xticklabel',t);
      if ( figure_row == 5 || i == size(cmb_total,1) )
         xlabel('Time (ms)');
      end
      set(gca,'ytick', ftick, 'yticklabel',f./1000);
      title(sprintf('2: %.0f-%.0f %.0fum %.2f Hz', chan2, model2, position2, w02));


      if ( mod(figure_row, 5) )
         figure_row = figure_row + 1;
         set(gcf,'position',[650 360 440 600]);
      else
         set(gcf,'position',[650 360 440 600]);
         figure_row = 1;
         figure;
      end

      cmap = brewmaps('rdbu', 21);
      colormap(cmap);

end

print_mfilename(mfilename);

return;




%    title(sprintf('%.0f:%.0f-%.0f %.0fum n_0%.0f w_0%.2f', ...
%       i, chan, model, position, n0, w0));
% 
%       orient tall;
%       suptitle(sprintf('%s site%.0f %s', exp, site, stim));
%       print_mfilename(mfilename);
% 
% return;









function [frequnique, ampunique, tc] = freq_resp_area_matrix(data, win)

freq = data.freq;
amp = data.amp;

frequnique = unique(freq);
ampunique = unique(amp);

tc = zeros(length(ampunique), length(frequnique));

for i = 1:length(freq)

   resp = data.resp{i};
   spetindex = find( resp>=win(1) & resp<=win(2) );

   if ( ~isempty(spetindex) )

      numspikes = length(spetindex);

      inda = find( amp(i) == ampunique); 
      indf = find( freq(i) == frequnique);

      tc(inda, indf) = tc(inda, indf) + numspikes;

   end

end % (for i)

return














