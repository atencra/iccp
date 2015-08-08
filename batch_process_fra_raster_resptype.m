function batch_process_fra_raster_resptype
% batch_process_fra_raster_resptype Temporal response type for FRAs
% 
%    batch_process_fra_raster_resptype searches through the current
%    directory for *-fracmb-pairs.mat files. For each file, the fra
%    variable is processed, and the response type is calculated. The
%    response type is estimated from the PSTH and raster, and provides a
%    measure of how phasic or tonic the response is.
%
%    Thus, this function is a wrapper to the function that does the
%    calculations, fra_raster_response_type.m

respfiles = dir( sprintf('*-fracmb-pairs.mat') );

for i = 1:length(respfiles)

   infile = respfiles(i).name;
   fprintf('Processing %s\n', infile);

   index = findstr(infile, '.mat');
   base_name = infile(1:index-1); % for later output file names

   outfile = sprintf('%s-resptype.mat', base_name); % output file name
   doutfile = dir(outfile); % see if it already exists

   if ( isempty(doutfile) ) % It doesn't exist, so process the data ...
      load(infile, 'fra');
      rt = fra_raster_response_type(fra);
      save(outfile, 'rt');
      fprintf('Saved data in %s\n\n', outfile);
   else % data already exists / already processed
      fprintf('Data exists in %s\n\n', outfile);
   end % if ( length(d)==0 )

end % (for i)

return;







