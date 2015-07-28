function iccp_batch_singleunit_spk2locator
% batch_singleunit_spk2locator Spike times to locator variable for single units
% 
%    Searches through the current directory and gets spk variables from
%    '*-strfcmb-pairs.mat' files. spk variables are struct array that
%    contain the spikes times for the neuron. One neuron's data is
%    contained in an element of an spk variable.
%
%    The spike times are then converted to binned spike trains using the
%    ripple stimulus triggers. This allows the binned spike times to be
%    aligned to the stimulus.
%
%    The output data is saved in '*-strfcmb-pairs-locator.mat' files.

% Need to first combine spike models before running this function ...

stimfolder = 'F:\20040217\stimuli_20040217';

dsfile = 'dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8-matrix.mat';
specfile = 'dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8.spr';
paramfile = 'dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8_param.mat';


pfiles = dir('*-strfcmb-pairs.mat');

for j = 1:length(pfiles)

   infile = pfiles(j).name;

   if ( 1 ) %length(spkfile) == 1 && length(trigfile) == 1 )

      fprintf('Getting Locator: SU for %s\n', infile);

      index = findstr(infile, '.mat');
      outfile = sprintf('%s-locator.mat',infile(1:index-1));
      fprintf('Outfile = %s\n', outfile);

      d = dir(outfile);

      if ( length(d)==0 ) % We haven't previously processed the data ...

         load(infile, 'spk');
         load(infile, 'trigger');
         spkloc = spk;
         exp = spkloc(1).exp;
         fs = spkloc(1).fs;

         fullspecfile = fullfile(stimfolder, specfile);
         fulldsfile = fullfile(stimfolder, dsfile);

         for k = 1:length(spkloc)

            sp = spkloc(k).spiketimes; % spiketimes are in ms
            spet = round(sp / 1000 * fs); % convert ms to sample number
            spet = spet(:)';
            trigger = trigger(:)';

            [taxis, faxis, locator, ~, ~] = ...
            get_locator_for_mid_analysis(fullspecfile, 0.2, 0.05, spet, trigger, fs);

            %             % The following will let you check the process to make sure
            %             % it's working.
            %             load(fulldsfile);
            %             [stim, resp] = get_mne_stim_resp(stimulus, locator, size(stimulus,1), 20, 1);
            %             sta = stim' * resp;
            %             figure; imagesc(reshape(sta,size(stimulus,1),20));

            spkloc(k).specfile = specfile;
            spkloc(k).paramfile = paramfile;
            spkloc(k).dsfile = dsfile;
            spkloc(k).taxis = taxis;
            spkloc(k).faxis = faxis;
            spkloc(k).locator = locator;
         end % (for k)

         fprintf('Saving Locator in %s\n\n', outfile);
         save(outfile, 'spkloc');
         clear('spk', 'trigger', 'spkloc');

      else
          fprintf('Outfile = %s already exists.\n', outfile);
      end % (if length(d)==0)

   end % (if)

end % (for j)

return;



