function spknew = iccp_isi_filtered_spk(spk, delay)
% ICCP_ISI_FIlTERED_SPK Remove ISI from spk spiketimes
%
% spknew = iccp_isi_filtered_spk(spk, delay)
%
% spk : struct array holding spiketimes. One element per neuron. Obtained
% from saveSpikeSort.m
%
% delay : minimum delay required between successive spikes. Spikes than
% occur less than delay ms apart are removed.
%
% spknew : new spk struct array with corrected spike times.


spknew = spk;

for j = 1:length(spk)

   spiketimes = spk(j).spiketimes;

   badindex = [];

   for i = 1:(length(spiketimes)-1)

      if ( abs( spiketimes(i) - spiketimes(i+1) ) < delay )
         badindex = [badindex i+1];
      end

   end % (for i)

   totalindex = 1:length(spiketimes);
   goodindex = setdiff(totalindex,badindex);
   newspiketimes = spiketimes(goodindex);

   spknew(j).spiketimes = newspiketimes;

end % (for j)











