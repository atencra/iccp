function iccp_batch_sta_fio_fitParams
% iccp_batch_sta_fio_fitParams Fit theoretical curve to nonlinearities in directory
%
% Gets all *-sta-fio.mat files in the current directory, calls 
% get_icc_fio_curve_fit_params on the fio data to obtain curve fit parameters.
%
% The curve fit has two parameters, a threshold and a transition variable.
% These paremeters describe the threshold to spiking for the neuron and
% then transition noise of the spike generating mechanism. If the neuron is
% a hard rectifier, then noise ~ 0, while values > 0 indicate less fidelity
% in coding auditory stimuli.
%
% Nomenclature: sta : spike triggered average
%               fio : function input/output - the nonlinearity
%               fitParams : parameters for curve fit
%
% Curve fit data are saved in files ending in *-sta_fio_fitparams.mat



d = dir('*-sta-fio.mat')
    
for n = 1:length(d)
    filename = d(n).name;
    load(filename, 'fio');
    
    index = findstr(filename, '-sta-fio');
    prefix = filename(1:index-1);
    newfile = sprintf('%s-sta-fio-fitparams.mat', prefix);

    fioFit = get_icc_fio_curve_fit_params(fio);
    
    save(newfile,'fioFit');
       
end



end
