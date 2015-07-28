function iccp_batch_strf_rtf
% iccp_batch_strf_rtf Estimate modulaton functions for STRFs in directory
% 
%    Processes all '*-pairs-strf-params.mat' files in the current directory
%    and then calculates ripple transfer function parameters. RTFs are the
%    Fourier transform of the STRF, and describe the temporal and spectral
%    modulation parameters that drive a neuron to fire.
%
%    Data are saved in '*-pairs-rtf-params.mat' files.

pfiles = dir('*-pairs-strf-params.mat');

if ( nargin == 0 )
   override = 0;
end

for j = 1:length(pfiles)

   infile = pfiles(j).name

   fprintf('Getting STRF Params for %s\n', infile);

   index = findstr(infile, '-strf-params.mat');
   outfile = sprintf('%s-rtf-params.mat',infile(1:index-1));
   fprintf('Outfile = %s\n', outfile);

   d = dir(outfile);

   if ( length(d)==0 ) % We haven't previously processed the data ...

      load(infile, 'params');

      rtf_params = strf_ripple_transfer_function(params);

      fprintf('Saving STRF Params in %s\n\n', outfile);
      save(outfile, 'rtf_params');
      clear('params', 'rtf_params');

   else
       fprintf('Outfile = %s already exists.\n', outfile);
   end % (if length(d)==0)

end % (for j)



return;



