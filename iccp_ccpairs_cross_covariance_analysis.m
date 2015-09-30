function ccpairs_new = iccp_ccpairs_cross_covariance_analysis(ccpairs)
% iccp_ccpairs_cross_covariance_analysis  Analyze peaks in cross-cov functions
%
%    ccpairs_new = iccp_ccpairs_cross_covariance_analysis(ccpairs)
%
%    Reanalyzes cross-covariance functions in the struct array ccpairs so 
%    that positive and negative peaks are evaluated. Why? Because the 
%    positive peaks represent neural excitation, while the negative peaks 
%    represent suppression. It is likely, then, that the characteristics 
%    of each type of peak are different.
%
%    ccpairs_new contains the processed data. This may appear a little
%    ad-hoc, mainly because this was an analysis that was thought of after
%    all initial processing had already taken place.


library('spike_train_cross_correlation');
library('stimbox');

for i = 1:length(ccpairs)

   fprintf('Processing #%.0f of %.0f\n', i, length(ccpairs) );

   dt = ccpairs(i).dt; % in ms

   spiketimes1 = ccpairs(i).spiketimes1;
   spiketimes2 = ccpairs(i).spiketimes2;

   delay = ccpairs(i).delay;
   r12 = ccpairs(i).r12;
   q12 = ccpairs(i).q12;
   conf_limit = ccpairs(i).conf_limit;

   exp = ccpairs(i).exp;
   [~, nsec] = stim_experiment_specfile(exp);
   duration = nsec * 1000; % duration of stimulus in ms

    [pd_pos, hw_pos, hw_pos_left, hw_pos_right, sigfeature_pos] = ...
        cc_crosscov_pd_hw(delay, q12, conf_limit);

   ccc_pos = cross_corr_coef(delay, r12, pd_pos, duration, dt, spiketimes1, spiketimes2);


   [pd_neg, hw_neg, hw_neg_left, hw_neg_right, sigfeature_neg] = ...
      cc_crosscov_pd_hw(delay, -q12, conf_limit);

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















