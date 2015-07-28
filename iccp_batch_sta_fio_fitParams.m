function iccp_batch_sta_fio_fitParams
% iccp_batch_sta_fio_fitParams Nonlinearity curve fits to directory data
%
% Finds all the *-sta-fio.mat files in the current directory, loads the
% nonlinearities contained in the file, then fits a parametric curve to
% each nonlinearity. The fit and parameters are then saved in
%
% *-sta-fio-fitparams.mat files.
%
% The theoretical curve is discussed in detail in Ringach and Malone
% (2009).


d = dir('*-sta-fio.mat');
    
for n = 1:length(d)
    filename = d(n).name;
    load(filename, 'fio');
    
    index = findstr(filename, '-sta-fio');
    prefix = filename(1:index-1);
    newfile = sprintf('%s-sta-fio-fitparams.mat', prefix);

    fioFit = iccp_fio_curve_fit_params(fio);
    save(newfile,'fioFit');
end



return;





