function vs_batch_strf_parameters(override)
% batch_singleunit_spk2locator Spike times to locator variable for single units
% 
%    batch_raw2thresh is run outside of the data folders holding all the .raw
%    files. It expects that within each folder there are .raw files for each
%    type of stimulus. It then processes each .raw file, finds multi-unit
%    data, and saves the data. It then moves on to th next folder.
% 
%    The stimuli that this function looks for have specific labels in the file
%    name. The stimuli must be ripple stimuli. Other stimuli will cause the
%    program to crash.
%
%    For experiments 2013-9-18, 2013-10-17, 2013-11-25:
%        stimuli =  {'rn11',  'rn41',  'rn12',  'rn42'};
% 
%    For experiment 2013-12-19:
%        stimuli =  {'rn1',  'rn4'};
%
%    If you have different stimuli, then you can run batch_raw2thresh(stimuli),
%    where stimuli is a cell array of strings holding different stimulus types.
%
%    This function expects there to be 16 channels of data. If there are
%    more or less then the function will need to be modified.
% 
%    Craig Atencio 11/20/13

% Need to first combine spike models before running this function ...


% something's not working right with this function

if ( nargin == 0 )
   override = 0;
end

pfiles = dir('*-strfcmb-pairs.mat');

for j = 1:length(pfiles)

   infile = pfiles(j).name;

   fprintf('Getting STRF Params for %s\n', infile);

   index = findstr(infile, '-strfcmb-pairs.mat');
   outfile = sprintf('%s-strfcmb-pairs-params.mat',infile(1:index-1));
   fprintf('Outfile = %s\n', outfile);

   d = dir(outfile);

   if ( length(d)==0 || override ) % We haven't previously processed the data ...

      load(infile, 'strf');
      load(infile, 'trigger');

      params = vs_strf_parameters(strf, trigger);

      fprintf('Saving STRF Params in %s\n\n', outfile);
      save(outfile, 'params');
      clear('strf', 'trigger', 'params');

   else
       fprintf('Outfile = %s already exists.\n', outfile);
   end % (if length(d)==0)

end % (for j)



return;



