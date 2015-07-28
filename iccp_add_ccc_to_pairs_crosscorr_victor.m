function ccpairs_new = iccp_add_ccc_to_pairs_crosscorr_victor(ccpairs)
% iccp_add_ccc_to_pairs_crosscorr_victor Add cross-correlation coef to cross-corr data struct
% 
% ccpairs_new = iccp_add_ccc_to_pairs_crosscorr_victor(ccpairs)
% 
% ccpairs : cross-correlation struct array. Each element corresponds to one pair of neurons.
% 
% ccpairs_new : same as ccpairs, but with the ccc added.
% 
% ccpairs_new = ccpairs;

for i = 1:length(ccpairs)

   exp = ccpairs(i).exp;

   nsec = 11*60;
   duration = nsec * 1000; % duration of stimulus in ms
   
%    fprintf('%s  nsec = %.0f  dur = %.0f\n', exp, nsec, duration);
%    pause
   
   delay = ccpairs(i).delay;
   r12 = ccpairs(i).r12;
   peakdelay = ccpairs(i).peakdelay;
   dt = ccpairs(i).dt;
   s1 = ccpairs(i).spiketimes1;
   s2 = ccpairs(i).spiketimes2;
   rho = ccpairs(i).rho;

   maxtime = max( [max(s1) max(s2)] );
   k = ceil( duration / dt ); % number of time bins in spike train

   index = find( peakdelay == delay );
   r12max = r12(index);

   if ( length(s2) > length(s1) )
      n1 = length(s1);
      n2 = length(s2);
   else
      n1 = length(s2);
      n2 = length(s1);
   end


   % Cross-Correlation Coefficient
   ccc = ( r12max - n1 * n2 / k ) / ( n1 - n1 * n2 / k );

   ccpairs_new(i).ccc = ccc;
   
end


return;


