function ccpairs_new = iccpairs_ccpairs_cross_covariance_analysis(ccpairs)
% iccpairs_ccpairs_cross_covariance_analysis  Metrics from cross-correlations
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


for i = 1:length(ccpairs)

   fprintf('Processing #%.0f of %.0f\n', i, length(ccpairs) );

   fsspk = ccpairs(i).fs;
   fsd = ccpairs(i).fsd; % Hz = 0.5 ms resolution
   dt = ccpairs(i).dt; % in ms

   spiketimes1 = ccpairs(i).spiketimes1;
   spiketimes2 = ccpairs(i).spiketimes2;
   n1 = ccpairs(i).n1;
   n2 = ccpairs(i).n2;

   delay = ccpairs(i).delay;
   r12 = ccpairs(i).r12;
   q12 = ccpairs(i).q12;
   conf_limit = ccpairs(i).conf_limit;

   [train1, train2] = corr_spiketimes2train(spiketimes1, spiketimes2, fsspk, fsd);

   maxdelay = 50;
   [delay, rab] = corr_crosscorr(train1, train2, maxdelay, dt);

   ca = corr_asi(delay, r12);
   cent = corr_centroid(delay, r12);


   exp = ccpairs(i).exp;
   [specfile, nsec] = get_specfile_for_exp(exp);
   duration = nsec * 1000; % duration of stimulus in ms

   [q12, conf_limit, stimdur] = get_cross_covariance(r12, n1, n2, dt, duration);


   [pd_pos, hw_pos, hw_pos_left, hw_pos_right, sigfeature_pos] = ...
      crosscov_params(delay, q12, conf_limit);

   ccc_pos = cross_corr_coef(delay, r12, pd_pos, duration, dt, spiketimes1, spiketimes2);


   [pd_neg, hw_neg, hw_neg_left, hw_neg_right, sigfeature_neg] = ...
      crosscov_params(delay, -q12, conf_limit);

   ccc_neg = cross_corr_coef_neg(delay, r12, -q12, pd_neg, duration, dt, spiketimes1, spiketimes2);

   % figure;
   % bar(delay, q12, 'k');
   % ccc_pos
   % ccc_neg
   % pause

   ccpairs(i).pd_pos = pd_pos;
   ccpairs(i).hw_pos = hw_pos;
   ccpairs(i).hw_pos_left = hw_pos_left;
   ccpairs(i).hw_pos_right = hw_pos_right;
   ccpairs(i).sigfeature_pos = sigfeature_pos;
   ccpairs(i).ccc_pos = ccc_pos;

   ccpairs(i).pd_neg = pd_neg;
   ccpairs(i).hw_neg = hw_neg;
   ccpairs(i).hw_neg_left = hw_neg_left;
   ccpairs(i).hw_neg_right = hw_neg_right;
   ccpairs(i).sigfeature_neg = sigfeature_neg;
   ccpairs(i).ccc_neg = ccc_neg;


% fprintf('\nPeak delay = %.2f\n', pd);
% fprintf('\nCentroid = %.2f\n', cent);
% fprintf('\nAsymmetry = %.2f\n', ca);
% fprintf('\nHalf-width = %.2f\n', hw);
% fprintf('\nSignificant? = %.2f\n', sigfeature);
% pause;

end % (for i)

ccpairs_new = ccpairs;


return;



function [pd, hw, hw_left_delay, hw_right_delay, sigfeature] = crosscov_params(delay, qab, conf_limit)
% crosscov_params Peak delay and half-width from cross-covariance function
% 
%    [pd, hw, hw_left_delay, hw_right_delay, sigfeature] = crosscov_params(delay, qab, conf_limit)
% 
%    delay : cross covariance delay
%    qab : cross covariance function
%    conf_limit : threshold for statistical significance
% 
%    pd : peak delay
%    hw : half width
%    hw_left_delay : left boundary of half width
%    hw_right_delay : right boundary of half width
%    sigfeature : 1 if a significant peak was present; 0 if not. A significant
%    peak is present if the peak and one cross-covariance bin to the left
%    or right of the peak both exceed the confidence limit. This is the
%    standard requirement that two consecutive bins meet statistical
%    significance.
% 
%    When significance is not reached, hw values are NaN.


qab(qab<0) = 0; % set negative values to 0

d = diff(delay);
d = d(1); % assuming equal bin widths
maxdelay = 20;
window = round(maxdelay / d);
ind0 = find(delay == 0);
indneg = find(delay<=- window);
indneg = max(indneg);
indpos = find(delay>=window);
indpos = min(indpos);

indcenter = indneg:indpos;
delaycenter = delay(indcenter);
qabcenter = qab(indcenter);

% Delay at peak of cross-covariance function
[temp, indmax] = max(qabcenter);
indpd = max(indmax);
pd = delaycenter(indpd);
peakvalue = qabcenter(indpd);


if ( pd > -maxdelay && pd < maxdelay )

   % process if there are two consecutive significant cross-covariance bins
   if ( peakvalue >= conf_limit && ...
      ( (qabcenter(indmax-1) >= conf_limit) || ...
      (qabcenter(indmax+1) >= conf_limit) )  ) 

      indpeakvalue = find(qab == peakvalue);
      delaypeakvalue = abs( delay(indpeakvalue) );
      index_delay = find(delaypeakvalue == min(min(delaypeakvalue)) );
      indpeakvalue = indpeakvalue(index_delay);

      index_leftmost = max([indpeakvalue-window 1]);
      index_left = index_leftmost:indpeakvalue;
      left = qab(index_left);
      delayleft = delay(index_left);
      indleft = find(left < 0.5 * peakvalue);

      index_rightmost = min([indpeakvalue+window length(qab)]);
      index_right = indpeakvalue:index_rightmost;
      right = qab(index_right);
      delayright = delay(index_right);
      indright = find(right < 0.5 * peakvalue);

      if ( ~isempty(indleft) && ~isempty(indright) )
         hw_left_delay = delayleft( max(indleft) );
         hw_right_delay = delayright( min(indright) );

      elseif ( ~isempty(indleft) && isempty(indright) )
         hw_left_delay = delayleft( max(indleft) );
         hw_right_delay = max(delayright);

      elseif ( isempty(indleft) && ~isempty(indright) )
         hw_left_delay = min(delayleft);
         hw_right_delay = delayright( min(indright) );

      else
         hw_left_delay = min(delayleft);
         hw_right_delay = max(delayright);

      end

      hw = hw_right_delay - hw_left_delay;
      sigfeature = 1;

   else
      sigfeature = 0;
      hw = nan;
      hw_left_delay = nan;
      hw_right_delay = nan;
   end


else
   sigfeature = 0;
   hw = nan;
   hw_left_delay = nan;
   hw_right_delay = nan;
end



% clf;
% hold on;
% ymax = max([max(qab) conf_limit]);
% range = ymax;
% bar(delay, qab, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
% plot([0 0], [0 ymax+0.1*range], 'k-');
% plot([min(delay) max(delay)], [0 0], 'k-');
% plot([min(delay) max(delay)], [conf_limit conf_limit], 'r-');
% plot([min(delay) max(delay)], 0.5*[peakvalue peakvalue], 'b--');
% tickpref;
% 
% plot([pd pd], [ylim], 'r-');
% plot([hw_left_delay hw_left_delay], [ylim], 'r-');
% plot([hw_right_delay hw_right_delay], [ylim], 'r-');
% title(sprintf('Significant = %.0f', sigfeature));
% pause(0.5);

return;






function ca = corr_asi(delay, rab)

d = diff(delay);
d = d(1); % assuming equal bin widths

window = round(20 / d);

ind0 = find(delay == 0);
indneg = find(delay<=- window);
indneg = max(indneg);
indpos = find(delay>=window);
indpos = min(indpos);

indcenter = indneg:indpos;

rabcenter = rab(indcenter);
delaycenter = delay(indcenter);


% Correlation Asymmetry Index from cross-correlation function
right = rab(ind0+1:indpos); 
left = rab(indneg:ind0-1); 
ca = ( sum(right) - sum(left) ) / ( sum(right) + sum(left) );

return;





function cent = corr_centroid(delay, rab)

d = diff(delay);
d = d(1); % assuming equal bin widths

window = round(20 / d);

ind0 = find(delay == 0);
indneg = find(delay<=- window);
indneg = max(indneg);
indpos = find(delay>=window);
indpos = min(indpos);

indcenter = indneg:indpos;

rabcenter = rab(indcenter);
delaycenter = delay(indcenter);

% Centroid of central portion of cross-correlation function
cent = sum( rabcenter .* delaycenter ) / sum(rabcenter);

return;









