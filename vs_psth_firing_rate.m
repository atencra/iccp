function [p, rt, edges] = vs_psth_firing_rate(trialspikes,dt,T)

% plot_psth Generates raster plot, PSTH, and stimulus spectrogram for a
% given response file and dt value

% [p, rt] = plot_psth(resp,dt)

% Inputs:
% ------------------------------
% 1) trialspikes: cell array. Each element contains the responses to one
%    presentation of a stimulus segment. The length of trialspikes is the
%    number of repetitions of the segment
% 2) dt: time bin size used to calculate the psth, in seconds.
% 3) T : total duration of the stimulus segment, in seconds. Default = 12.
%
% Outputs:
% ------------------------------
% 1) p: psth binned at resolution dt
% 2) rt: firing rate from p, in spikes/s

if ( nargin == 2 )
   T = 12;
end

edges = 0:dt:T;

all_spiketimes = [];

for j = 1:length(trialspikes)
   spikes = trialspikes{j}./1000; % spike times in ms; put them in seconds
   all_spiketimes = [all_spiketimes spikes];
end

p = histc(all_spiketimes, edges);

rt = p ./ length(trialspikes) ./ dt;
     
return



