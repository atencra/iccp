function [idt, dt] = vs_info_time_bin(trialspikes, dt, T)
% NEEDS TO BE MODULARIZED

% info_calculation: Calculates information values and nmse values from an
% STRF response file. (saved as resp_final)

% [info_nmse, sigma_nmse, nmse] = info_calculation(resp, graphs, popgraphs)

% Inputs:
% -----------------------------
% 1) resp: STRF response file, saved as -*final-resp. 
% 2) graphs: boolean input. 1 if graphs are desired, 0 if not.
% 3) popgraphs: boolean input. 1 if population graphs are desired, 0 if
%       not.

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

%   IF POPGRAPHS is 1 (plots are for all the neurons)
%     1) Plots the population histogram of normalized mean square error, 
%             information, and picked dt values.

%  Finally, performs pairwise analysis on information, sigma values, and
%  generates scatterplots of those.

idt = zeros(1,length(dt));

for n = 1:length(dt)

   [p, rt] = vs_psth_firing_rate(trialspikes, dt(n), T);
   rbar = mean(rt);

  index = find( rt>0 );
  idt(n) = dt(n) * 1/T .* sum( (rt(index)./rbar).*log2( rt(index)./rbar) );

end

return;





