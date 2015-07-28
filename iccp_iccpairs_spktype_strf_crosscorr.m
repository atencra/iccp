function vs_iccpairs_spktype_strf_crosscorr(override)
% vs_iccpairs_spktype_strf_crosscorr Correlation between pairs
%
% vs_iccpairs_spktype_strf_crosscorr(override)
% -----------------------------------------------------------
%
% Reads in data from pairs analysis and estimates the cross-correlation
% between the spike trains of each neuron in a pair. Data is read from 
% files in a directory that have already been processed. The data files
% are:
%
% *-strfcmb-pairs.mat, 
% *-strfcmb-pairs-params.mat,
% *-strfcmb-pairs-ptparams.mat, 
% *-strfcmb-pairs-rtf-params.mat,
% *-strfcmb-pairs-sta-fio.mat
%
% Each file contains a struct array. The number of struct elements in each
% struct array for each file must match.
%
% If the elements match, the cross-correlation between the spike trains for
% each pair are estimated, and the results are saved in a struct array that
% contains all the data in the individual files as well as the new
% cross-correlation data.
%
% The new struct array is then saved in a file with the ending: 
% *-strfcmb-pairs-ccpairs.mat
%
% If the override input argument is not supplied, then the function first
% checks to make sure that there is not already an output *-ccpairs.mat
% file. If there is, then the function skips processing the data for that
% set of files. If override = 1, then the function runs and writes over any
% *-ccpairs.mat file that is present.
%
% caa 8/15/14


narginchk(0,1);
if ( nargin == 0 )
   override = 0;
end

pfiles = dir('*-strfcmb-pairs.mat');

for j = 1:length(pfiles)

   strffile = pfiles(j).name;

   fprintf('Getting STRF Params for %s\n', strffile);

   index = findstr(strffile, '.mat');
   basename = strffile(1:index-1);
   outfile = sprintf('%s-ccpairs.mat', basename );
   fprintf('Outfile = %s\n', outfile);

   d = dir(outfile);

   if ( length(d)==0 || override ) % We haven't previously processed the data ...

      load(strffile, 'strf', 'spk', 'trigger');

      paramsfile = sprintf('%s-params.mat', basename );
      load(paramsfile, 'params');

      ptparamsfile = sprintf('%s-ptparams.mat', basename );
      load(ptparamsfile, 'ptparams');

      rtfparamsfile = sprintf('%s-rtf-params.mat', basename );
      load(rtfparamsfile, 'rtf_params');

      rtfparamsfile = sprintf('%s-sta-fio.mat', basename );
      load(rtfparamsfile, 'fio');

      rtfparamsfile = sprintf('%s-sta-fio-fitparams.mat', basename );
      load(rtfparamsfile, 'fioFit');




      % check to make sure that each struct has the same number of elements
      nFiles = length(strf);
      totalFiles = length(strf) + length(ptparams) + length(params) + ...
          length(rtf_params) + length(fio);

      if ( 5*nFiles ~= totalFiles )
         error('Mismatch between number of files for different analyses.');
      end

      [rfpairs] = ...
         get_pairs_spktype_spiketrains(spk, strf, trigger, ptparams, params, rtf_params, fio);

      [ccpairs] = get_pairs_spktype_fsu_rsu_crosscorr(rfpairs);

      fprintf('Saving CCPairs in %s\n\n', outfile);
      save(outfile, 'rfpairs', 'ccpairs');
      clear('spk', 'strf', 'trigger', 'params', 'ptparams', 'rtf_params', 'fio');

   else
       fprintf('Outfile = %s already exists.\n', outfile);
   end % (if length(d)==0)

end % (for j)

return;



%------------------- Function Definitions ------------------------

function [datastr] = get_pairs_spktype_spiketrains(spk, strf, trigger, ptparams, params, rtf_params, fio)
% For pairs of neurons from the same channel, get the data, and save it
% in a struct array for later processing
%
% The data that we'll save are the experiment, site, channel, model,
% depth, position, raw spike trains, and wave form durations,


p = 0.01;

datastr = [];

chan = [spk.chan];
chan_unique = unique(chan);

for i = 1:length(chan_unique)

   index = find( chan_unique(i) == chan );

   if ( index > 1 )

      cmb = nchoosek(index, 2); % determine all possible pairwise combinations

      [nr, nc] = size(cmb);

      for j = 1:nr

      index1 = cmb(j,1);
      index2 = cmb(j,2);

      % First get the spike shape data

      % Neuron #1, phase 1 and 2
      %                begphase1 = spktype(index1).begphase;
      %                endphase1 = spktype(index1).endphase;

      % Neuron #2, phase 1 and 2
      %                begphase2 = spktype(index2).begphase;
      %                endphase2 = spktype(index2).endphase;
      %                
      %                dur1 = begphase1 + endphase1;
      %                dur2 = begphase2 + endphase2;


      % assign the data to a temporary struct variable
      datatemp.exp = spk(index1).exp;
      datatemp.site = spk(index1).site;
      datatemp.chan = spk(index1).chan;

      datatemp.model1 = spk(index1).model;
      datatemp.model2 = spk(index2).model;

      datatemp.depth = spk(index1).depth;
      datatemp.position = spk(index1).position;
      datatemp.stim = spk(index1).stim;
      datatemp.atten = spk(index1).atten;

      datatemp.header1 = spk(index1).header;
      datatemp.header2 = spk(index2).header;

      datatemp.waveform1 = spk(index1).waveform;
      datatemp.waveform2 = spk(index2).waveform;

      datatemp.spiketimes1 = spk(index1).spiketimes;
      datatemp.spiketimes2 = spk(index2).spiketimes;

      datatemp.outliers1 = spk(index1).outliers;
      datatemp.outliers2 = spk(index2).outliers;

      datatemp.meanwaveform1 = spk(index1).meanwaveform;
      datatemp.meanwaveform2 = spk(index2).meanwaveform;

      datatemp.fs = spk(index1).fs;

      %                datatemp.dur1 = dur1;
      %                datatemp.dur2 = dur2;


      % Get the pure tone parameters
      cf1 = ptparams(index1).cf;
      cf2 = ptparams(index2).cf;

      q1 = ptparams(index1).q;
      q2 = ptparams(index2).q;

      latency1 = ptparams(index1).latency;
      latency2 = ptparams(index2).latency;

      offset = 0;		
      latency1 = latency1 + offset;
      latency2 = latency2 + offset;

      datatemp.cf1 = cf1;
      datatemp.cf2 = cf2;

      datatemp.q1 = q1;
      datatemp.q2 = q2;

      datatemp.latency1 = latency1;
      datatemp.latency2 = latency2;


      % Now get the STRF similarity data
      rf1 = strf(index1).rfcontra;
      rf2 = strf(index2).rfcontra;

      n01 = strf(index1).n0contra;
      n02 = strf(index2).n0contra;

      mdb = strf(index1).mdb;

      fs = strf(index1).fs;

      stimdur = ( trigger(end) - trigger(1) ) / fs;

      soundtype = strf(index1).stim;

      rfsig1 = significant_strf(rf1, p, n01, mdb, stimdur, soundtype);
      rfsig2 = significant_strf(rf2, p, n02, mdb, stimdur, soundtype);

      [cc, std1, std2] = similarity_index(rfsig1, rfsig2);

      datatemp.strfsimilarity = cc;


      % Get the firing rate, phase locking index, and separability index
      fr1 = params(index1).w0;
      fr2 = params(index2).w0;

      pli1 = params(index1).pli;
      pli2 = params(index2).pli;

      eigvals1 = params( index1 ).eigvals;
      si1 = eigvals1(1) / (sum(eigvals1)+eps);

      eigvals2 = params( index2 ).eigvals;
      si2 = eigvals2(1) / (sum(eigvals2)+eps);


      datatemp.fr1 = fr1;
      datatemp.fr2 = fr2;

      datatemp.pli1 = pli1;
      datatemp.pli2 = pli2;

      datatemp.sepindex1 = si1;
      datatemp.sepindex2 = si2;



      % Get modulation parameters
      btmf1 = rtf_params(index1).best_tmf{end};
      btmf2 = rtf_params(index2).best_tmf{end};

      bsmf1 = rtf_params(index1).best_xmf{end};
      bsmf2 = rtf_params(index2).best_xmf{end};

      twc3db1 = rtf_params(index1).twc3db{end}(end);
      twc3db2 = rtf_params(index2).twc3db{end}(end);

      swc3db1 = rtf_params(index1).swc3db{end}(end);
      swc3db2 = rtf_params(index2).swc3db{end}(end);


      datatemp.btmf1 = btmf1;
      datatemp.btmf2 = btmf2;

      datatemp.bsmf1 = bsmf1;
      datatemp.bsmf2 = bsmf2;

      datatemp.twc3db1 = twc3db1;
      datatemp.twc3db2 = twc3db2;

      datatemp.swc3db1 = swc3db1;
      datatemp.swc3db2 = swc3db2;


      % Get the nonlinearity parameters
      fio1 = fio(index1);
      fio2 = fio(index2);

      [sta1Mn, pspkx1Mn, x1Mn, y1max] = get_pairs_sta_fio_stats(fio1);
      [sta2Mn, pspkx2Mn, x2Mn, y2max] = get_pairs_sta_fio_stats(fio2);

      [x12, i1, i2] = intersect(x1Mn, x2Mn);

      pspkx1 = pspkx1Mn(i1);
      pspkx2 = pspkx2Mn(i2);

      r = corrcoef(pspkx1, pspkx2); % correlation between nonlinearities
      fiosimilarity = r(1,2);

      asi1 = fio_asi(x12, pspkx1);
      asi2 = fio_asi(x12, pspkx2);

      datatemp.fiosimilarity = fiosimilarity;
      datatemp.fioasi1 = asi1;
      datatemp.fioasi2 = asi2;


      % assign the temporary struct to a larger struct array

      datastr = [datastr datatemp];
      clear('datatemp');

      end % (for j)

   end % (if)

end % (for i)

return;




function [fsufsu, fsursu, rsursu] = get_pairs_spktype_fsu_rsu_spiketrains(datastr)

% get the spike waveform durations
dur1 = [datastr.dur1];
dur2 = [datastr.dur2];

% get indices for cell type pairs
[index_fsu_fsu, index_rsu_rsu, index_fsu_rsu, index_rsu_norsu, ...
    index_fsu_norsu, index_norsu_norsu] = fsu_rsu_indices(dur1, dur2);

% separate data into three main types of pairs
fsufsu = datastr(index_fsu_fsu);
rsursu = datastr(index_rsu_rsu);
fsursu = datastr(index_fsu_rsu);

return;




function [ccstr] = get_pairs_spktype_fsu_rsu_crosscorr(unitstr)

ccstr = unitstr;

for i = 1:length(ccstr)

   fprintf('Processing #%.0f of %.0f\n', i, length(ccstr) );

   fsspk = ccstr(i).fs;

   spiketimes1 = ccstr(i).spiketimes1;
   spiketimes2 = ccstr(i).spiketimes2;
   maxspiketime = max([ max(spiketimes1) max(spiketimes2) ] );
   spiketimes1 = [spiketimes1 maxspiketime];
   spiketimes2 = [spiketimes2 maxspiketime];

   spet1 = spiketimes1 ./ 1000 .* fsspk;
   spet2 = spiketimes2 ./ 1000 .* fsspk;

   fsd = 2000; % Hz = 0.5 ms resolution
   dt = 1 / fsd * 1000; % in ms

   train1 = spet2train(spet1, fsspk, fsd);
   train2 = spet2train(spet2, fsspk, fsd);


   % Compute cab = E[xa(t+u)xb(t)]
   [r12, delay] = xcorr(train1, train2, 100); % make max lag = 50 * 0.5 ms = 25 ms
   r12 = real(r12);
   r12(find(r12 < 0)) = 0;
   delay = delay .* dt; % delay in ms, not bin number
   n1 = length(spet1);
   n2 = length(spet2);

   % save the data to the struct array
   ccstr(i).fsd = fsd;
   ccstr(i).dt = dt;
   ccstr(i).n1 = n1; % the number of spikes
   ccstr(i).n2 = n2; % the number of spikes
   ccstr(i).delay = delay;
   ccstr(i).r12 = r12;

   [q12, conf_limit, stimdur] = get_cross_covariance(r12, n1, n2, dt);
   [pd, cent, ca, hw, sigfeature] = ...
      get_pairs_cross_covariance_params(delay, q12, r12, conf_limit);

   rho = corrstrength(delay, r12, n1, n2, stimdur); % correlation coefficient

%    delay = ccpairs(i).delay;
%    r12 = ccpairs(i).r12;
%    peakdelay = ccpairs(i).peakdelay;
%    dt = ccpairs(i).dt;
%    s1 = ccpairs(i).spiketimes1;
%    s2 = ccpairs(i).spiketimes2;
%    rho = ccpairs(i).rho;
%    maxtime = max( [max(spiketimes1) max(spiketimes2)] );

   k = length(train1); %ceil( duration / dt ); % number of time bins in spike train
   index = find( pd == delay );
   r12max = r12(index);

   if ( length(spiketimes2) > length(spiketimes1) )
      n1 = length(spiketimes1);
      n2 = length(spiketimes2);
   else
      n1 = length(spiketimes2);
      n2 = length(spiketimes1);
   end

   % Cross-Correlation Coefficient
   ccc = ( r12max - n1 * n2 / k ) / ( n1 - n1 * n2 / k );

   ccstr(i).q12 = q12;
   ccstr(i).conf_limit = conf_limit;
   ccstr(i).rho = rho;
   ccstr(i).peakdelay = pd;
   ccstr(i).centroid = cent;
   ccstr(i).asymmetry = ca;
   ccstr(i).halfwidth = hw;
   ccstr(i).significant = sigfeature;
   ccstr(i).ccc = ccc;
    
   clf;
   subplot(2,1,1);
   bar(delay, r12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
   tickpref;

   subplot(2,1,2);
   upper95qab = conf_limit;
   lower95qab = -conf_limit;
   hold on;
   ymin = min([min(q12) lower95qab]);
   ymax = max([max(q12) upper95qab]);
   range = ymax-ymin;
   bar(delay, q12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
   plot([0 0], [ymin-0.1*range ymax+0.1*range], 'k-');
   plot([min(delay) max(delay)], [0 0], 'k-');
   plot([min(delay) max(delay)], [upper95qab upper95qab], 'r-');
   plot([min(delay) max(delay)], [lower95qab lower95qab], 'r-');
   tickpref;

   % fprintf('\nPeak delay = %.2f\n', pd);
   % fprintf('\nCentroid = %.2f\n', cent);
   % fprintf('\nAsymmetry = %.2f\n', ca);
   % fprintf('\nHalf-width = %.2f\n', hw);
   % fprintf('\nSignificant? = %.2f\n', sigfeature);

   pause(0.1);
    
end % (for i)


return;




