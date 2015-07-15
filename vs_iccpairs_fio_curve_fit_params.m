function fioFit = vs_iccpairs_fio_curve_fit_params(fio)
% get_icc_fio_curve_fit_params  Curve fit to midbrain nonlinearity data
%
% fioFit = get_icc_fio_curve_fit_params(fio) fits a theoretical curve
% to the nonlinearity data in fio. The nonlinearity data that is fit
% is for the STA.
%
% caa 5/24/12



fields = {'locator', 'itrain'};
fio = rmfield(fio,fields);
fioFit = fio;

for i = 1:length(fio)

   x = fio(i).xbins;
   bins = x;
   pspkx = fio(i).pspkx;
   fx = double(nanmean([pspkx{1} pspkx{2} pspkx{3} pspkx{4}], 2));
   pspkx = fx; % reassigned for later use in the function

   pspk = fio(i).pspk;
   f0 = double(nanmean([pspk{1} pspk{2} pspk{3} pspk{4}], 2));
   pspk = f0;


   % remove NaN's from the nonlinearity
   index = ~isnan(fx);
   x = x(index);
   bins = bins(index);
   fx = fx(index);
   pspkx = pspkx(index);


   dt = 0.001;

%    fx = pspkx ./ dt;
%    f0 = pspk ./ dt;
%    fx = fx - f0;
%    fx(fx<0) = 0;


   % Don't use the last point for curve fitting since it is subject 
   % to estimation error
   x = x(1:end-1);
   fx = fx(1:end-1);

   % Fit function from Ringach and Malone (2007)
   a0 = [1 1 1]; % initial guess for the params
   [fitParams, resNorm] = lsqcurvefit(@ringach_malone_func, a0, x, fx);


   % The curve fit at higher resolution
   xFit = linspace(min(x), max(x), 100);
   fxFit = ringach_malone_func(fitParams, xFit);


   % Get goodness of fit metrics for the curve fit
   fxTemp = ringach_malone_func(fitParams, x);
   nmse = gfit(fx, fxTemp, '2'); % gets goodness of fit measure, nmse
   r2 = gfit(fx, fxTemp, '8'); % gets goodness of fit measure, coef of determination


   % Calculate nonlinearity asymmetry index (ASI)
   right = pspkx(bins>0.1);
   left = pspkx(bins<-0.1);
   asi = ( sum(right) - sum(left) ) / ( sum(right) + sum(left) );   


   % Calculate the peak rate and skewness of the nonlinearity
   peakRate = 1./ dt .* max(pspkx);
   skew = skewness(pspkx);


   % Calculate the projection value when the nonlinearity crosses
   % the average background rate, pspk
   index = find(bins > -1);
   bins_right = bins(index);
   pspkx_right = pspkx(index);
   bins_right = bins_right(1:end-2);
   pspkx_right = pspkx_right(1:end-2);
   xi = linspace(min(bins_right), max(bins_right), 1000); % projection values
   yi = interp1(bins_right, pspkx_right, xi, 'spline'); % interpolated nonlinearity

   index = find(yi > pspk, 1);
   threshold = xi(index);    % Calculate the crossing point


   % Save the data to a struct; a '0' after the field name means that
   % the parameter was calculated for the STA
   fioFit(i).dt = dt;
   fioFit(i).x_sta = x;
   fioFit(i).fx_sta = fx;
   fioFit(i).xFit_sta = xFit;
   fioFit(i).fxFit_sta = fxFit;
   fioFit(i).fitParams_sta = fitParams;
   fioFit(i).nmse_sta = nmse;
   fioFit(i).r2_sta = r2;
   fioFit(i).asi_sta = asi;
   fioFit(i).peakRate_sta = peakRate;
   fioFit(i).skew_sta = skew;
   fioFit(i).threshold_sta = threshold;


%    clf
%    hold on;
%    plot(x, fx, 'ko', 'markerfacecolor', 'k');
%    plot(xFit, fxFit, 'r-');
%    plot([fitParams(2) fitParams(2)], [0 max(fx)], 'k-');
%    xmax = max([max(x) abs(min(x))]);
%    xmax = xmax + 2*xmax*0.05;
%    xlim([-xmax xmax]);
%    ylimit = get(gca,'ylim');
%    set(gca,'ylim', [-0.1*max(fx) max(ylimit)]);
%    set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    legend('Data', 'Fit', 'location', 'northwest');
%    xlabel('Similarity between STA and Stimulus (SD)');
%    ylabel('Firing Rate (Hz)');
%    title(sprintf('#%.0f of %.0f: NMSE = %.3f, R2 = %.3f', ...
%       i, length(fio), nmse, r2));
% 
% pause



end % (for i)


%fioFit


return;



