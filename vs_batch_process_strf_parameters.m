function vs_batch_process_strf_parameters
% vs_batch_process_strf_parameters Calculate STRF parameters from files in a directory
% 
%    batch_process_acute_rat_tm(stimfolder, stimuli) goes through all the directories
%    in an experimental folder and processes the data. This function should be
%    run inside a folder that contains individual folders for each penetration.
%    Inside the penetration folders are *-thresh.mat and *-trig.mat files, which
%    will be used to process the data.
% 
%    Craig Atencio 
%    8/13/14


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







