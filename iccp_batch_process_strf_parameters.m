function iccp_batch_process_strf_parameters
% iccp_batch_process_strf_parameters Calculate STRF parameters from files in a directory
% 
%    Goes through all the '*-strfcmb-pairs.mat' files in the current
%    directory and calculates STRF parameters, such separability index,
%    response precision, ripple transfer function, etc.
% 
%    Data are save in '*-strfcmb-pairs-params.mat' files


files = dir( sprintf('*-strfcmb-pairs.mat') );

for i = 1:length(files)

   infile = files(i).name;
   fprintf('Processing %s\n', infile);
   index = findstr(infile, '.mat');
   base_name = infile(1:index-1); % for later output file names
   outfile = sprintf('%s-params.mat', base_name); % output file name
   doutfile = dir(outfile); % see if it already exists

   if ( isempty(doutfile) ) % It doesn't exist, so process the data ...
      load(infile, 'strf', 'trigger');
      params = vs_strf_parameters(strf,trigger);
      save(outfile, 'params');
      fprintf('Saved data in %s\n\n', outfile);
   else % data already exists / already processed
      fprintf('Data exists in %s\n\n', outfile);
   end % if ( length(d)==0 )

end % (for i)

return;







