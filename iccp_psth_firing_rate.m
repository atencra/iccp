function [p, rt, edges] = iccp_psth_firing_rate(trialspikes,dt,T)
% iccp_psth_firing_rate PSTH, firing rate from spike times
% 
%     [p, rt, edges] = iccp_psth_firing_rate(trialspikes,dt,T)
% 
%     trialspikes : cell array of spike times. length(trialspikes) = # stimulus
%     repetitions.
% 
%     dt : bin size of PSTH, in seconds.
%     T : duration of PSTH, in seconds.
%     p : PSTH
% 
%     rt : time varying firing rate
%     edges : PSTH bin edges
% 
%     Output values are used as inputs to information calculations.



if ( nargin == 2 )
   T = 12; % duration of ICC repeated stimulus
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



