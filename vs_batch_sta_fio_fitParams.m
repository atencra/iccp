function vs_batch_sta_fio_fitParams
% get_sta_fio_fitParams Calculates nonlinearity curve fit parameters

% Gets all *-sta-fio.mat files, calls get_icc_fio_curve_fit_params on the
% fio data to obtain curve fit parameters.

% Saves output data as files ending in *-sta_fio_fitparams.mat


d = dir('*-sta-fio.mat');
    
for n = 1:length(d)
    filename = d(n).name;
    load(filename, 'fio');
    
    index = findstr(filename, '-sta-fio');
    prefix = filename(1:index-1);
    newfile = sprintf('%s-sta-fio-fitparams.mat', prefix);

    fioFit = vs_iccpairs_fio_curve_fit_params(fio);
    save(newfile,'fioFit');
end



return;





