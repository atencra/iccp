function iccp_batch_add_strf_tm_field_value
% iccp_batch_add_strf_tm_field_value Add temporal modulation value STRF struct array
% 
%     For all '*-strfcmb-pairs.mat', loads the STRF struct array and adds the TM field.
%     It give the field a value of 500, which is the maximum stimulus temporal
%     modulation frequency.
% 
%     This was missing from the data structures and is needed for later analyses.


pfiles = dir('*-strfcmb-pairs.mat');

for j = 1:length(pfiles)

   infile = pfiles(j).name;

   fprintf('Changing STRF TM field value for %s\n', infile);

   load(infile, 'strf');
   load(infile, 'spk');
   load(infile, 'trigger');

   for i = 1:length(strf)
      strf(i).tm = 500;
   end % (for i)

   fprintf('Saving STRF in %s\n\n', infile);
   save(infile, 'strf', 'spk', 'trigger');
   clear('strf', 'spk', 'trigger');

end % (for j)


return;



