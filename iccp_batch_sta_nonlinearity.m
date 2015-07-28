function iccp_batch_sta_nonlinearity(stimfolder)
% iccp_batch_sta_nonlinearity Estimate nonlinearity for files in directory
 %
% iccp_batch_sta_nonlinearity(stimfolder)
% --------------------------------------------------------
%
% This function goes through the current directory, opens *-strf-pairs.mat
% files, and computes a binned spike train that is a downsampled version of
% the original spike train. The bin size of the spike train is usually ~5
% ms. The stimulus is then obtained, and the nonlinearity for each
% receptive field is calcualted. Nonlinearity data is saved in 
%
% stimfolder is the absolute path to the directory holding the stimuli that
% were used to create the earlier non-downsampled STRFs. This folder will
% hold a downsampled version of the original stimulus, which will be used
% to estimate the nonlinearity.
%
% Data are saved in *-sta-fio.mat files.
%

if ( nargin ~= 1 )
   error('You need to specify the stimulus folder containing the .spr files.');
end

d = dir('*-pairs-locator.mat');

t1 = 0;
t2 = 0.1;

for n = 1:length(d)

   infile = d(n).name;
   fprintf('Getting locator for %s\n', infile);
   index = findstr(infile, '.mat');
   outfile = sprintf('%s-sta-fio.mat',infile(1:index-1));

   fprintf('Outfile = %s\n', outfile);
   d = dir(outfile);

   if ( length(d)==0 ) % We haven't previously processed the data ...

       load(infile);

       fprintf('Processing %s\n', filename);

       fio = spkloc; % spkloc will hold the locator for later analysis
       exp = spkloc(1).exp;
       [specfile, paramfile, dsfile] = get_downsampled_specfile(exp, 'dmr');

       specfile = spkloc(1).specfile;
       specfile = fullfile(stimfolder, specfile); % total path to spr file

       stim_matrix_file = fullfile(stimfolder, dsfile);
       load(stim_matrix_file);

       for i = 1:length(spkloc)

          locator = spkloc(i).locator;

          [sta, itrain] = get_sta_train_from_locator(locator, stimulus, 30);

          [nf, ntrials] = size(stimulus);

          index_freq = 1:nf; % use entire frequency axis for calculations

          [x0train, x0test, x0train_locator, x0test_locator] = ...
             train_test_projection(sta, locator, stimulus, index_freq);

          [x0bins, pspk, px0, px0spk, pspkx0] = ...
             proj_prob_dist(x0train, x0train_locator); % for the nonlinearity

          fio(i).sta = sta;
          fio(i).itrain = itrain;
          fio(i).xbins = x0bins;
          fio(i).pspk = pspk;
          fio(i).px = px0;
          fio(i).pxspk = px0spk;
          fio(i).pspkx = pspkx0;

       end % (for i)

       % Save spkloc to a new file if it doesn't already exist

       if ( isempty( dir(outfile) ) )
          fprintf('Outfile %s\n', outfile);
    %       eval(['save ' outfile ' fio']);
          save(outfile, 'fio');
       else
          fprintf('\nThe file %s already exists.\n', outfile);
       end

       clear spkloc fio stimulus* locator spkloc

    end

end % (for n)












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








function get_pairs_sta_nonlinearity_folder(pathname)

d = dir('*-pairs-locator.mat');

t1 = 0;
t2 = 0.1;

for n = 1:length(d)

   filename = d(n).name;
   load(filename);

   index = findstr(filename, 'locator.mat');
   outfile = [filename(1:index-1) 'sta-fio.mat'];

   fprintf('Processing %s\n', filename);

   fio = spkloc; % spkloc will hold the locator for later analysis

%   spkloc = rmfield( spkloc, {'header' 'waveform' 'spiketimes' ...
%      'outliers' 'fs' 'meanwaveform'} ); % remove fields we won't need

   exp = spkloc(1).exp;

   specfile = spkloc(1).specfile;

   specfile = [pathname '\' specfile]; % total path to spr file

   index = findstr(specfile, '.spr');
   stim_matrix_file = [specfile(1:index-1) '-matrix.mat'];

   load(stim_matrix_file);

   for i = 1:length(spkloc)

      locator = spkloc(i).locator;

      [sta, itrain] = get_sta_train_from_locator(locator, stimulus, 30);

      [nf, ntrials] = size(stimulus);

      index_freq = 1:nf; % use entire frequency axis for calculations

      [x0train, x0test, x0train_locator, x0test_locator] = ...
         train_test_projection(sta, locator, stimulus, index_freq);

      [x0bins, pspk, px0, px0spk, pspkx0] = ...
         proj_prob_dist(x0train, x0train_locator); % for the nonlinearity

      fio(i).sta = sta;
      fio(i).itrain = itrain;
      fio(i).xbins = x0bins;
      fio(i).pspk = pspk;
      fio(i).px = px0;
      fio(i).pxspk = px0spk;
      fio(i).pspkx = pspkx0;

   end % (for i)

   % Save spkloc to a new file if it doesn't already exist

   if ( isempty( dir(outfile) ) )
      fprintf('Outfile %s\n', outfile);
%       eval(['save ' outfile ' fio']);
      save(outfile, 'fio');
   else
      fprintf('\nThe file %s already exists.\n', outfile);
   end

   clear spkloc fio stimulus* locator spkloc

end % (for n)


return;











