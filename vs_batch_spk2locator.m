function vs_batch_spk2locator(stimfolder)
% get_pairs_locator - locator for pairs of neurons
% 
% get_pairs_locator(stimfolder)
% --------------------------------------------------------
%
% This function goes through the current directory, opens *-strf-pairs.mat
% files, and computes a binned spike train that is a downsampled version of
% the original spike train. The bin size of the spike train is usually ~5
% ms, allowing for later MID and nonlinearity analysis.
%
% stimfolder is the absolute path to the directory holding the stimuli that
% were used to create the earlier non-downsampled STRFs.
%
% caa 4/12/11

% F:\20031119\stimuli_20031119
% F:\20031124\stimuli_20031124
% F:\20040217\stimuli_20040217

if ( nargin ~= 1 )
   error('You need to specify the stimulus folder containing the .spr files.');
end

% Get the files with the data
d = dir('*-strfcmb-pairs.mat');


% These values are set, but not used
t1 = 0;
t2 = 0.1;

for n = 1:length(d)

   filename = d(n).name;
   load(filename);

   index = findstr(filename, '.mat');
   outfile = [filename(1:index-1) '-locator.mat'];

   fprintf('\nProcessing %s\n', filename);
   fprintf('Outfile = %s\n', outfile);

   spkloc = spk; % spkloc will hold the locator for later analysis

   spkloc = rmfield( spkloc, {'header' 'waveform' 'spiketimes' ...
      'outliers' 'fs' 'meanwaveform'} ); % remove fields we won't need

   exp = spk(1).exp;
   stim = spk(1).stim;

   specfile = get_pairs_downsampled_specfile(exp, stim);

% specfile
% pause

   specfile_path = fullfile(stimfolder,specfile); % total path to spr file

   for i = 1:length(spk)

      spiketimes = spk(i).spiketimes; % spiketimes in ms
      fsspk = spk(i).fs;
      spet = spiketimes / 1000 * fsspk; % spiketimes in sample number

      [taxis, faxis, locator, numspikes, averate, locator_ipsi] = ...
         get_locator_for_mid_analysis(specfile_path, t1, t2, ...
            spet, trigger, fsspk); 

      spkloc(i).taxis = taxis;
      spkloc(i).faxis = faxis;
      spkloc(i).locator = locator;
      spkloc(i).specfile = specfile;

   end % (for i)

   % Save spkloc to a new file if it doesn't already exist

   if ( isempty( dir(outfile) ) )
      save(outfile, 'spkloc');
   else
      fprintf('\nThe file %s already exists.\n', outfile);
   end

   clear spk spkloc trigger

end % (for n)


return;







