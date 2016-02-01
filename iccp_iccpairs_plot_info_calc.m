function iccp_iccpairs_plot_info_calc(iresp)
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

% Things to plot:
% 1. Raster
% 2. PSTH
% 3. Information vs. time bin size
% 4. Information vs. data set size for multiple estimates
% 5. SD vs. data set size


nrows = 4;
ncols = 4;

chan = [iresp.chan];
chan_unique = unique(chan);
chan_unique = 16;

close all;

for i = 1:length(chan_unique)

   index = find( chan_unique(i) == chan );

   if ( length(index)>1 )

      cmb = nchoosek(index, 2); % determine all possible pairwise combinations

      [nr, nc] = size(cmb);

      for j = 1:nr

         index1 = cmb(j,1);
         index2 = cmb(j,2);

         iresp1 = iresp(index1);
         iresp2 = iresp(index2);

         dt1 = iresp1.dt;
         dt2 = iresp2.dt;

         idt1 = iresp1.idt;
         idt2 = iresp2.idt;

         dtoptim1 = iresp1.dtoptim;
         dtoptim2 = iresp2.dtoptim;

         [dtoptim1_temp] = vs_info_pick_dt(idt1, dt1);
         [dtoptim2_temp] = vs_info_pick_dt(idt2, dt2);



         % Make pairwise plots
         figure;

         % Raster plots
         subplot(nrows,ncols,1);
         vs_iccpairs_plot_iresp_raster(iresp1);
         ylabel('Repetition #');

         subplot(nrows,ncols,2);
         vs_iccpairs_plot_iresp_raster(iresp2);


         % PSTH plots
         subplot(nrows,ncols,5);
         vs_iccpairs_plot_iresp_psth(iresp1);
         ylabel('Firing Rate (sp/s)');
         xlabel('Time (s)');

         subplot(nrows,ncols,6);
         vs_iccpairs_plot_iresp_psth(iresp2);
         xlabel('Time (s)');


         % Info vs dt
         subplot(nrows,ncols,3);
         vs_plot_info_pick_dt(idt1, dt1);
         ylabel('Information (bits)');
         xlabel('Bin Size (s)');

         subplot(nrows,ncols,4);
         vs_plot_info_pick_dt(idt2, dt2);
         xlabel('Bin Size (s)');



         % Info vs. dt plots
         subplot(nrows,ncols,7);
         vs_iccpairs_plot_iresp_dt(iresp1);
         ylabel('Information (bits)');
         xlabel('1/Duration (1/s)');

         subplot(nrows,ncols,8);
         vs_iccpairs_plot_iresp_dt(iresp2);
         xlabel('1/Duration (1/s)');


         % Mean info vs. dt plots
         subplot(nrows,ncols,11);
         vs_iccpairs_plot_iresp_dt_infomn(iresp1);
         ylabel('Information (bits)');
         xlabel('1/Duration (1/s)');

         subplot(nrows,ncols,12);
         vs_iccpairs_plot_iresp_dt_infomn(iresp2);
         xlabel('1/Duration (1/s)');


         % SD info vs. dt plots
         subplot(nrows,ncols,15);
         vs_iccpairs_plot_iresp_dt_infosd(iresp1);
         ylabel('Information SD (bits)');
         xlabel('1/sqrt(Duration)');
         

         subplot(nrows,ncols,16);
         vs_iccpairs_plot_iresp_dt_infosd(iresp2);
         xlabel('1/sqrt(Duration)');
         
         set(gcf,'position', [680    30   730 640]);
         print_mfilename(mfilename);

% pause




         if ( length(chan_unique) == 1 )

            % Make pairwise plots
            figure;

            % Raster plots
            subplot(6,2,1);
            vs_iccpairs_plot_iresp_raster(iresp1);
            ylabel('Repetition #');

            subplot(6,2,2);
            vs_iccpairs_plot_iresp_raster(iresp2);


            % PSTH plots
            subplot(6,2,3);
            vs_iccpairs_plot_iresp_psth(iresp1);
            ylabel('Firing Rate (sp/s)');
            xlabel('Time (s)');

            subplot(6,2,4);
            vs_iccpairs_plot_iresp_psth(iresp2);
            xlabel('Time (s)');


            % Info vs dt
            subplot(6,2,5);
            vs_plot_info_pick_dt(idt1, dt1);
            ylabel('Information (bits)');
            xlabel('Bin Size (s)');

            subplot(6,2,6);
            vs_plot_info_pick_dt(idt2, dt2);
            xlabel('Bin Size (s)');



            % Info vs. dt plots
            subplot(6,2,7);
            vs_iccpairs_plot_iresp_dt(iresp1);
            ylabel('Information (bits)');
            xlabel('1/Duration (1/s)');

            subplot(6,2,8);
            vs_iccpairs_plot_iresp_dt(iresp2);
            xlabel('1/Duration (1/s)');


            % Mean info vs. dt plots
            subplot(6,2,9);
            vs_iccpairs_plot_iresp_dt_infomn(iresp1);
            ylabel('Information (bits)');
            xlabel('1/Duration (1/s)');

            subplot(6,2,10);
            vs_iccpairs_plot_iresp_dt_infomn(iresp2);
            xlabel('1/Duration (1/s)');


            % SD info vs. dt plots
            subplot(6,2,11);
            vs_iccpairs_plot_iresp_dt_infosd(iresp1);
            ylabel('Information SD (bits)');
            xlabel('1/sqrt(Duration)');


            subplot(6,2,12);
            vs_iccpairs_plot_iresp_dt_infosd(iresp2);
            xlabel('1/sqrt(Duration)');

            set(gcf,'position', [680 30 321 955]);
            print_mfilename(mfilename);

         end % (if length())


      end % (for j)

   end % (if/else)

 end % (for i)

return;





function vs_iccpairs_plot_iresp_raster(iresp)

raster = iresp.trialspikes;
total_dur = iresp.total_dur;
exp = iresp.exp;
site = iresp.site;
chan = iresp.chan;
model = iresp.model;

for i = 1:length(raster)
   hold on;
   plot(raster{i}./1000, i * ones(size(raster{i})), 'k.', 'markersize', 1, 'markerfacecolor', 'k');  
end % (for i)
xlim([0 total_dur]);
ylim([0 length(raster)+1]);

for j = 1:length(model)
   m(j) = num2str(model(j));
end % (for j)
title(sprintf('s%.0f c%.0f m%s', site, chan, m ));
set(gca,'ytick', [1 75 150], 'yticklabel', [1 75 150]);
tickpref;

return; 




function vs_iccpairs_plot_iresp_psth(iresp)

dt = 0.010;
total_dur = iresp.total_dur;
rbar = iresp.rbar;
hold on;
[p, rt, edges] = vs_psth_firing_rate(iresp.trialspikes, dt, total_dur);
hb = bar(edges, rt, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
plot(xlim, [rbar rbar], 'r-');
xlim([0 total_dur])
ylim([0 max(rt)]);
box off;
tickpref;

return;





function vs_iccpairs_plot_iresp_dt(iresp)

fracdur = iresp.frac_dur;
ifracdur = iresp.ifracdur;

ivals_total = [];
for ii = 1:length(ifracdur)
   ivals = ifracdur{ii};
   ivals_total = [ivals_total ivals(:)'];
   hold on;
   plot(ones(size(ivals)).*1./fracdur(ii), ivals, 'ko', 'markersize', 3);
end

set(gca,'xscale', 'log');
xlim([min(1./fracdur) max(1./fracdur)]);
ylim([0 1.1*max(ivals_total)]);
frac_tick = [1/12 1/6 1/3 1/1.201];
frac_label = {'1/12', '1/6', '1/3', '1/1.2'};
set(gca,'xtick', frac_tick, 'xticklabel', frac_label);
tickpref;

return;




function vs_iccpairs_plot_iresp_dt_infomn(iresp)

fracdur = iresp.frac_dur;
ifracdur = iresp.ifracdur;
ifracdurmn = iresp.ifracdurmn;
ifracdursd = iresp.ifracdursd;

% frac_tick = 1 ./ fracdur;
% frac_tick = fliplr(frac_tick(:)');
% for i = 1:length(fracdur)
%    frac_label{i} = sprintf('1/%.2f',fracdur(i));
% end % (for i)
% temp = cell(size(frac_label));
% for i = 1:length(frac_label)
%    temp{i} = frac_label{length(frac_label)-i+1};
% end
% frac_label = temp;
% index = [1 9 12 13 length(frac_tick)];
% frac_tick = frac_tick(index);
% frac_label = frac_label(index);

frac_tick = [1/12 1/6 1/3 1/1.201];
frac_label = {'1/12', '1/6', '1/3', '1/1.2'};

% fracdur
% 1./fracdur
% frac_label
% temp

plot(1./fracdur, ifracdurmn, 'ko', 'markersize', 3);
set(gca,'xscale', 'log');
xlim([min(1./fracdur) max(1./fracdur)]);
ylim([0 1.25*max(ifracdurmn)]);
set(gca,'xtick', frac_tick, 'xticklabel', frac_label);
tickpref;
box off;

return;





function vs_iccpairs_plot_iresp_dt_infosd(iresp)

fracdur = iresp.frac_dur;
ifracdur = iresp.ifracdur;
ifracdurmn = iresp.ifracdurmn;
ifracdursd = iresp.ifracdursd;

x = 1./sqrt(fracdur);
y = ifracdursd;
beta = polyfit(x, y, 1);
xfit = linspace(min(x), max(x), 500);
yfit = polyval(beta, xfit);
hold on;
plot(x, y, 'ko', 'markersize', 3);
plot(xfit, yfit, 'k-');
xlim([sqrt(min(1./fracdur)) sqrt(max(1./fracdur))]);
ylim([0 1.1*max(ifracdursd)]);
tickpref;

return;






function [dtoptim] = vs_plot_info_pick_dt(idt, dt)

dt_error = [];
idt_error = [];

for i = 1:length(dt)-25

   index = [i:length(dt)];
   dtsubset = dt(index);
   idtsubset = idt(index);

   beta = polyfit(dtsubset, idtsubset, 1);
   idtsubset_fit = polyval(beta, dtsubset);
   idtfit = polyval(beta, dt);

   e = abs( (idtsubset - idtsubset_fit) ./ idtsubset );
   idt_error = [idt_error max(e)];
   dt_error = [dt_error min(dtsubset)];

   fitstr(i).dt = dt;
   fitstr(i).idt = idt;
   fitstr(i).index = index;
   fitstr(i).dtsubset = dtsubset;
   fitstr(i).idtsubset = idtsubset;
   fitstr(i).beta = beta;
   fitstr(i).idtsubset_fit = idtsubset_fit;
   fitstr(i).idtfit = idtfit;

end % (for i )

index = find(idt_error > 0.2);
index = max(index);
dtoptim = dt_error(index);

hold on;
plot(dt, idt, 'ko', 'markersize', 3);
% plot(fitstr(index).dtsubset, fitstr(index).idtsubset_fit, 'k-', 'linewidth', 2);
plot(fitstr(index).dt, fitstr(index).idtfit, 'r-', 'linewidth', 2);
plot([dtoptim dtoptim], ylim, 'r-', 'linewidth', 2);
xlim([0-0.001 max(dt)]);
tickpref;

return;







