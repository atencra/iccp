function plot_all_nonlinearities
% plot_all_nonlinearities: Function that plots out all nonlinearity
% according to the curve fit parameters stored in *-sta-fio-fitparams.mat.
% Iterates over all *-sta-fio-fitparams.mat files and plots out the curves
% that they describe.

d = dir('*-sta-fio-fitparams.mat');
dL = length(d);
% maxnmse = [];

for n = 1:dL
    filename = d(n).name
    load(filename, 'fioFit');
    
    L = length(fioFit);
    %     nmse_values = [];
    
    for i = 1:L
        
        x = fioFit(i).x_sta;
        fx = fioFit(i).fx_sta;
        xFit = fioFit(i).xFit_sta;
        fxFit = fioFit(i).fxFit_sta;
        fitParams = fioFit(i).fitParams_sta;
        nmse = fioFit(i).nmse_sta;
        r2 = fioFit(i).r2_sta;
        
        if nmse < 0.6055    %maximum nmse value
            figure
            hold on;
            plot(x, fx, 'ko', 'markerfacecolor', 'k');
            plot(xFit, fxFit, 'r-');
            plot([fitParams(2) fitParams(2)], [0 max(fx)], 'k-');
            xmax = max([max(x) abs(min(x))]);
            xmax = xmax + 2*xmax*0.05;
            xlim([-xmax xmax]);
            ylimit = get(gca,'ylim');
            set(gca,'ylim', [-0.1*max(fx) max(ylimit)]);
            set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
            legend('Data', 'Fit', 'location', 'northwest');
            xlabel('Similarity between STA and Stimulus (SD)');
            ylabel('Firing Rate (Hz)');
            title(sprintf('#%.0f of %.0f: NMSE = %.3f, R2 = %.3f, Theta = %.2f, Sigma = %.2f',...
                i, L, nmse, r2, fitParams(2), fitParams(3)));
            
            pause;
            close
            %         nmse_values = [nmse_values;nmse];
        end
        
    end
    %     maxnmse = [maxnmse ; max(nmse_values)];
    
end
% maxnmse

end