function newspiketimes = iccp_isi_filtered_spiketimes(spiketimes, delay)
% iccp_filtered_spiketimes Remove ISI from spiketimes
%
% spknew = iccp_isi_filtered_spiketimes(spiketimes, delay)
%
% spiketimes : vector of spike times, in ms.
%
% delay : minimum delay required between successive spikes. Spikes than
% occur less than delay ms apart are removed.
%
% newspiketimes : spiketimes with no ISI < delay.


badindex = [];

for i = 1:(length(spiketimes)-1)

  if ( abs( spiketimes(i) - spiketimes(i+1) ) < delay )
     badindex = [badindex i+1];
  end

end % (for i)

totalindex = 1:length(spiketimes);
goodindex = setdiff(totalindex,badindex);
newspiketimes = spiketimes(goodindex);

return;






