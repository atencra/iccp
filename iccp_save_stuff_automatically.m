function save_stuff_automatically

d = dir('*-strfcmb-pairs.mat');

for n = 1:length(d)
    
    filename = d(n).name;
	s = load(filename, 'strf', 'trigger');
    strf = s.strf;
    trigger = s.trigger;
    
    index = findstr(filename, '-strfcmb-pairs');
    prefix = filename(1:index-1);
    newfile1 = sprintf('%s-pairs-strf-params.mat', prefix);
    
    params = strf_parameters(strf,trigger);
   
%     save(newfile1,'params') 
    
    newfile2 = sprintf('%s-pairs-rtf-params.mat', prefix);
    
    rtf_params = strf_ripple_transfer_function(params);
    
    save(newfile2,'rtf_params') 
end
   
