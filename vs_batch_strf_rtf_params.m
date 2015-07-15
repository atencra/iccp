function vs_batch_strf_rtf_params
% vs_batch_strf_rtf_params Calculate STRF parameters from files in a directory
% 
%    vs_batch_strf_rtf_params goes through all the strf parameter files
%    in the current directory and calculates ripple transfer function and 
%    modulation transfer functions. 
%
%    Craig Atencio 
%    8/13/14


files = dir( sprintf('*-strfcmb-pairs-params.mat') );

for i = 1:length(files)

   infile = files(i).name;
   fprintf('Processing %s\n', infile);
   index = findstr(infile, '-params.mat');
   base_name = infile(1:index-1); % for later output file names
   outfile = sprintf('%s-rtf-params.mat', base_name); % output file name
   doutfile = dir(outfile); % see if it already exists

   if ( isempty(doutfile) ) % It doesn't exist, so process the data ...
      load(infile, 'params');
      rtf_params = strf_ripple_transfer_function(params);
      save(outfile, 'rtf_params');
      fprintf('Saved data in %s\n\n', outfile);
   else % data already exists / already processed
      fprintf('Data exists in %s\n\n', outfile);
   end % if ( length(d)==0 )

end % (for i)

return;







