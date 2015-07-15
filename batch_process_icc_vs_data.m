function batch_process_icc_vs_data(stimfolder)
% batch_process_acute_rat_tm Calculate FRAs, TM responses, and STRFs from multiunit data.
% 
%    batch_process_acute_rat_tm(stimfolder, stimuli) goes through all the directories
%    in an experimental folder and processes the data. This function should be
%    run inside a folder that contains individual folders for each penetration.
%    Inside the penetration folders are *-thresh.mat and *-trig.mat files, which
%    will be used to process the data.
% 
%    stimfolder : absolute path to folder holding all the stimuli for the
%    experiments.
%       Example: stimfolder = 'F:\2013-9-18-trained-rat-tm20-acute\stimuli_2013-9-18';
% 
%    stimuli : cell array holding strings of stimulus types that were presented
%
%    For experiments 2013-9-18, 2013-10-17, 2013-11-25:
%        stimuli =  {'fra', 'tm1',  'rn11',  'rn41',  'tm2',   'rn12',  'rn42'};
% 
%    For experiment 2013-12-19:
%        stimuli =  {'fra1', 'fra2',  'tm', 'rn1',  'rn4'};
%
%    Note: you don't have to run all the stimuli at once. You can run
%    individual stimuli. In this case, you might set stimuli = {'fra1'},
%    for example.
%
%    Craig Atencio 
%    12/5/13

F:\20040217\stimuli_20040217

dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min.bin                          
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min.spr                          
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min.wav                          
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8-matrix.mat         
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8.spr                
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8_param.mat          
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_param.mat                    
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-12sec.bin                          
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-12sec.spr                          
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-12sec_150repeats_stim2_44khz.wav   
dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-12sec_param.mat                    



2003-11-19-site1-2366um-50db-dmr1-fs18115-locator.mat             
2003-11-19-site1-2366um-50db-dmr1-fs18115-pairs-rtf-params.mat    
2003-11-19-site1-2366um-50db-dmr1-fs18115-pairs-spktype.mat       
2003-11-19-site1-2366um-50db-dmr1-fs18115-pairs-strf-params.mat   
2003-11-19-site1-2366um-50db-dmr1-fs18115-sta-fio-fitparams.mat   
2003-11-19-site1-2366um-50db-dmr1-fs18115-sta-fio.mat             
2003-11-19-site1-2366um-50db-dmr1-fs18115-strf-ptparams.mat       
2003-11-19-site1-2366um-50db-dmr1-fs18115-strf.mat                
2003-11-19-site1-2366um-50db-dmr1-fs18115-strfcmb-pairs.mat   



respfiles = dir( sprintf('20*-*-site*-*um-*db-dmrrep-fs*-final-resp.mat') );

for i = 1:length(respfiles)

   infile = respfiles(i).name;
   fprintf('Processing %s\n', infile);
   index = findstr(infile, '.mat');
   base_name = infile(1:index-1); % for later output file names
   outfile = sprintf('%s-info.mat', base_name); % output file name
   doutfile = dir(outfile); % see if it already exists

   if ( isempty(doutfile) ) % It doesn't exist, so process the data ...
      load(infile, 'resp');
      iresp = vs_info_calc(resp);
      save(outfile, 'iresp');
      fprintf('Saved data in %s\n\n', outfile);
   else % data already exists / already processed
      fprintf('Data exists in %s\n\n', outfile);
   end % if ( length(d)==0 )

end % (for i)

return;







