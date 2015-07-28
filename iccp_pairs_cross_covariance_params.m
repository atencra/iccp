function [pd, cent, ca, hw, sigfeature] = vs_pairs_cross_covariance_params(...
    delay, qab, rab, conf_limit)
% get_pairs_cross_covariance_params  Metrics from cross-correlations
% 
%    [pd, cent, ca, hw, sigfeature] = ...
%          get_pairs_cross_covariance_params(delay, qab, rab, conf_limit)
% 
%    Inputs:
%    ----------------------------------------------
%    delay : lag of cross-covariance function, in ms
%    qab : cross-covariance function
%    rab : cross-correlation function
%    conf_limit : confidence limit. Used to assess significance
% 
%    Outputs:
%    ----------------------------------------------
%    pd : delay at which peak in cross-covariance function occurs
%    cent : centroid of cross-correlation function
%    ca : asymmetry of cross-correlation function
%    hw : half width of cross-covariance function
%    sigfeature : significance of function. Two consecutive cross-covariance 
%       function values must exceed conf_limit for significance
%
%


upper95qab = conf_limit;
lower95qab = -conf_limit;

d = diff(delay);
d = d(1); % assuming equal bin widths

window = round(20 / d);

ind0 = find(delay == 0);
indneg = find(delay<=- window);
indneg = max(indneg);
indpos = find(delay>=window);
indpos = min(indpos);

%[indneg ind0 indpos]

indcenter = indneg:indpos;
qabcenter = qab(indcenter);
rabcenter = rab(indcenter);
delaycenter = delay(indcenter);


% Correlation Asymmetry Index from cross-correlation function
right = rab(ind0+1:indpos); 
left = rab(indneg:ind0-1); 
ca = ( sum(right) - sum(left) ) / ( sum(right) + sum(left) );


% Centroid of central portion of cross-correlation function
cent = sum( rabcenter .* delaycenter ) / sum(rabcenter);

% Delay at peak of cross-covariance function
[temp, indmax] = max(qabcenter);
indpd = max(indmax);
pd = delaycenter(indpd);


% Make sure there are at least two consecutive significant 
% cross-covariance feature values
indqabsig = find( qabcenter >= upper95qab | qabcenter <= lower95qab );
delaysig = delaycenter(indqabsig);
diffdelaysig = diff(delaysig);

% only process if there are two consecutive significant cross-covariance
% features and if there are more than 10 coincident spikes in the peak
if ( ~isempty( find(diffdelaysig==d) ) && rabcenter(indpd) > 10 ) 
                                          
   peakvalue = qabcenter(indpd);
   indpeakvalue = find(qab == peakvalue);

   if ( length(indpeakvalue) > 1 ) % find peak nearest 0 delay
      indtemp = find(abs(indpeakvalue-ind0) == min(abs(indpeakvalue-ind0)) );
      indpeakvalue = indpeakvalue(indtemp);
   end

   if ( indpeakvalue < window+1 || indpeakvalue+window > length(qab) )
      sigfeature = 0;
      hw = nan;
   else
      left = qab(indpeakvalue-window:indpeakvalue);
      delayleft = delay(indpeakvalue-window:indpeakvalue);

      right = qab(indpeakvalue:indpeakvalue+window);
      delayright = delay(indpeakvalue:indpeakvalue+window);

      indleft = find(left < 0.5*peakvalue);
      indright = find(right < 0.5*peakvalue);

      if ( ~isempty(indleft) && ~isempty(indright) )
         hw = delayright(min(indright)) - delayleft(max(indleft));
      elseif ( ~isempty(indleft) && isempty(indright) )
         hw = max(delayright) - delayleft(max(indleft));
      elseif ( isempty(indleft) && ~isempty(indright) )
         hw = delayright(min(indright)) - min(delayleft);
      else
         hw = max(delayright) - min(delayleft);
      end
      sigfeature = 1;
   end

else
   sigfeature = 0;
   hw = nan;
end

return;

