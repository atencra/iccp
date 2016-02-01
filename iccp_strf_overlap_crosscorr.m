function [ccoverlap] = iccp_strf_overlap_crosscorr(ccpairs)
% iccp_strf_overlap_crosscorr STRF overlap for neuron pairs
%
% [ccoverlap] = iccp_strf_overlap_crosscorr
% -----------------------------------------------------------
%
% The function searches through '*-strfcmb-pairs-orig-isi-fixed.mat' and
% finds the neurons that were recorded from the 
% same channel, and then calculates the correlation between the STRFs
% for each pair of neurons.
%
% ccpairs : struct array holding correlation data/parameters. Each element
% of ccpairs holds the data for one pair of neurons.
%
% ccoverlap : struct array holding the overlap data for each pair of neurons.
% ccoverlap has a field, named "overlap", that contains the calculations.
% The other fields of ccoverlap are identical to those in ccpairs.



library('strfbox');
library('pairs');

narginchk(1,1);


% First find the cross-covariance functions that only have a positive peak:
[sigPosOnly, sigPosNeg, sigNegOnly] = iccp_ccpairs_sigfeature_index(ccpairs);



% Make reduced ccpairs struct from significant, positive only peak data:
ccpairsPos = ccpairs(sigPosOnly);



% Load find STRF data in files, estimate STRF overlap
[ccoverlap] = iccp_get_strf_overlap_from_files_ccpairs(ccpairsPos);


return;





%------------------- Function Definitions ------------------------

function [ccoverlap] = iccp_get_strf_overlap_from_files_ccpairs(ccpairs)
% iccp_get_strf_overlap_from_files_ccpairs Use files/ccpairs to get STRF overlap data
% 
%     [ccoverlap] = iccp_get_strf_overlap_from_files_ccpairs(ccpairs)
%     --------------------------------------------------------------
%     ccpairs : struct array holding cross-covariance data. Each element holds the
%     data for one pair of ICC neurons.
% 
%     ccoverlap : struct array holding the STRF overlap/correlation data for each
%     element in ccpairs.
% 
%     For each pair of neurons in ccpairs, the file name that should hold the 
%     STRF data is constructed. This file is then loaded from the current
%     directory, and the STRF overlap is estimated.
% 
%     This means that the function needs to be run inside the folder that holds
%     the STRF data files.


ccoverlap = ccpairs;

for ii = 1:length(ccpairs)

    exp = ccpairs(ii).exp;
    site = ccpairs(ii).site;
    depth = ccpairs(ii).depth;
    atten = ccpairs(ii).atten;
    stim = ccpairs(ii).stim;

    filePattern = sprintf('%s-site%.0f-%.0fum-%.0fdb-%s-fs*-strfcmb-pairs-isi-fixed.mat',...
                    exp, site, depth, atten, stim);

    strfFile = dir(filePattern);

    if length(strfFile) == 1
        load(strfFile.name);
        fprintf('%s\n', strfFile.name);


        % Get the STRFs for the ccpairs element
        [strf1, strf2] = pairs_get_strf_pair_for_crosscorr(strf, ccpairs(ii));

        if ( ~isempty(strf1) && ~isempty(strf2) )

            % Now get the STRF similarity data
            taxis = strf1.taxis;
            rf1 = strf1.rfcontra;
            rf2 = strf2.rfcontra;

            n01 = strf1.n0contra;
            n02 = strf2.n0contra;

            mdb = strf1.mdb;
            fs = strf1.fs;
            stimdur = ( trigger(end) - trigger(1) ) / fs;
            soundtype = strf1.stim;
            p = 0.01;

            rfsig1 = significant_strf(rf1, p, n01, mdb, stimdur, soundtype);
            rfsig2 = significant_strf(rf2, p, n02, mdb, stimdur, soundtype);

            [overlap] = strf_overlap(rfsig1, rfsig2, taxis);

            ccoverlap(ii).overlap = overlap;

        else
            ccoverlap(ii).overlap = [];
        end

    else
        fprintf('%s does not exist.\n', strfFile.name);
        ccoverlap(ii).overlap = [];

    end


    clear('spk', 'strf', 'trigger');

end % (for ii)
    

%     % construct file name from ccpairs(i) element.
% 
%     dstrf = dir('*-strfcmb-pairs-orig-isi-fixed.mat');
% 
% 
%     p = 0.01;
% 
%     datastr = [];
% 
%     for n = 1:length(dstrf)
% 
%         fprintf('Processing #%.0f of %.0f: %s\n', n, length(dstrf), dstrf(n).name);
% 
%         filename = dstrf(n).name;
%         s = load(filename, 'spk', 'strf', 'trigger');
%         spk = s.spk; 
%         strf = s.strf;
%         trigger = s.trigger;
%         chan = [spk.chan];
%         chan_unique = unique(chan);
% 
% 
%         index = findstr(filename, '-orig');
%         prefix = filename(1:index-1);
% 
%         dspktype = dir(sprintf('%s-spktype.mat',prefix));
%         filename = dspktype.name;
%         s = load(filename, 'spktype');
%         spktype = s.spktype;
% 
%         if ( length(spktype) ~= length(spk) )
%             error('spktype and strf are not the same length.');
%         end
% 
% 
%         for i = 1:length(chan_unique)
% 
%             index = find( chan_unique(i) == chan );
% 
%             cmb = nchoosek(index, 2); % determine all possible pairwise combinations
% 
%             [nr, nc] = size(cmb);
% 
%             for j = 1:nr
% 
%                 index1 = cmb(j,1);
%                 index2 = cmb(j,2);
% 
%                 % First get the spike shape data
% 
%                 % Neuron #1, phase 1 and 2
%                 begphase1 = spktype(index1).begphase;
%                 endphase1 = spktype(index1).endphase;
% 
%                 % Neuron #2, phase 1 and 2
%                 begphase2 = spktype(index2).begphase;
%                 endphase2 = spktype(index2).endphase;
% 
%                 dur1 = begphase1 + endphase1;
%                 dur2 = begphase2 + endphase2;
% 
% 
%                 % assign the data to a temporary struct variable
%                 datatemp.exp = spk(index1).exp;
%                 datatemp.site = spk(index1).site;
%                 datatemp.chan = spk(index1).chan;
% 
%                 datatemp.model1 = spk(index1).model;
%                 datatemp.model2 = spk(index2).model;
% 
%                 datatemp.depth = spk(index1).depth;
%                 datatemp.position = spk(index1).position;
%                 datatemp.stim = spk(index1).stim;
%                 datatemp.atten = spk(index1).atten;
% 
%                 datatemp.header1 = spk(index1).header;
%                 datatemp.header2 = spk(index2).header;
% 
%                 datatemp.waveform1 = spk(index1).waveform;
%                 datatemp.waveform2 = spk(index2).waveform;
% 
%                 datatemp.spiketimes1 = spk(index1).spiketimes;
%                 datatemp.spiketimes2 = spk(index2).spiketimes;
% 
%                 datatemp.outliers1 = spk(index1).outliers;
%                 datatemp.outliers2 = spk(index2).outliers;
% 
%                 datatemp.meanwaveform1 = spk(index1).meanwaveform;
%                 datatemp.meanwaveform2 = spk(index2).meanwaveform;
% 
%                 datatemp.fs = spk(index1).fs;
% 
%                 datatemp.dur1 = dur1;
%                 datatemp.dur2 = dur2;
% 
% 
% 
%                 % Now get the STRF similarity data
%                 taxis = strf(index1).taxis;
%                 rf1 = strf(index1).rfcontra;
%                 rf2 = strf(index2).rfcontra;
% 
%                 n01 = strf(index1).n0contra;
%                 n02 = strf(index2).n0contra;
% 
%                 mdb = strf(index1).mdb;
% 
%                 fs = strf(index1).fs;
% 
%                 stimdur = ( trigger(end) - trigger(1) ) / fs;
% 
%                 soundtype = strf(index1).stim;
% 
%                 rfsig1 = significant_strf(rf1, p, n01, mdb, stimdur, soundtype);
%                 rfsig2 = significant_strf(rf2, p, n02, mdb, stimdur, soundtype);
% 
%                 [cc, std1, std2] = similarity_index(rfsig1, rfsig2);
% 
%                 datatemp.strfsimilarity = cc;
% 
% 
%                 [overlap] = strf_overlap(rfsig1, rfsig2, taxis);
% 
%                 datatemp.overlap = overlap;
%                 datastr = [datastr datatemp];
%                 clear('datatemp');
% 
%             end % (for j)
% 
%         end % (for i)
% 
%     end % (for n)
% 
% end % (for i)

return;





function iccp_ccpairs_statistics(ccpairs)
% iccp_ccpairs_statistics Statistical results for Cross-cov, STRF, FIO parameters
% 
%    iccp_ccpairs_statistics(ccpairs) shows analysis results for the ICC
%    paired neuron analysis. The analysis assesses cross-covariance
%    functions, STRFs, and nonlinearities.
%
%    Data are grouped into Pos, Pos+Neg, and Neg categories, where Pos
%    indicates the cross-cov function had only a positive deflection,
%    Pos+Neg means there was a significant Pos and Neg deflection, and Neg
%    means there was only a Neg deflection.
%
%    There are no graphs displayed. This function is mainly to help in
%    gathering statistics for the paper.
% 
%    ccpairs : struct array, with each element holding the data for a pair
%    of ICC neurons.
% 


[sigPosOnly, sigPosNeg, sigNegOnly] = iccp_ccpairs_sigfeature_index(ccpairs);

iccp_ccpairs_crosscorr_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)

iccp_ccpairs_crosscorr_strf_si_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);

% iccp_sigfeature_to_fr_pli_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);
% 
% iccp_sigfeature_to_cf_q_lat_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);
% 
% iccp_sigfeature_to_bmf_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);
% 
% iccp_sigfeature_to_fio_si_asi_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);
% 
% iccp_sigfeature_to_fiofit_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly);


return;



function [sigPosOnly, sigPosNeg, sigNegOnly] = iccp_ccpairs_sigfeature_index(ccpairs)


pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);

sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);


if ( isfield(ccpairs,'pd_pos') && isfield(ccpairs,'pd_neg') )
    pdNeg = [ccpairs.pd_neg];
    pdNeg = abs(pdNeg);

    sigNeg = [ccpairs.sigfeature_neg];

    sigPosOnly = sigPos & ~sigNeg;

    sigNeg0 = logical(sigNeg & pdNeg > 0); % for positive/negative cross-cov
    sigNeg15 = logical(sigNeg & pdNeg > 1.5); % for negative only cross-cov

    sigNegOnly = ~sigPos & sigNeg15;
    sigPosNeg = sigPos & sigNeg0;

    fprintf('\n');
    fprintf('Positive only peaks: %.0f\n', sum(sigPosOnly) );
    fprintf('Negative only peaks: %.0f\n', sum(sigNegOnly) );
    fprintf('Positive and Negative peaks: %.0f\n', sum(sigPosNeg) );
    fprintf('Any peaks: %.0f\n', sum([sigPosOnly(:); sigNegOnly(:); sigPosNeg(:)]) );
    fprintf('\n');
end


if ( ~isfield(ccpairs,'pd_pos') && isfield(ccpairs,'peakdelay') )
    sigPosOnly = sigPos;
    sigNegOnly = false(size(sigPos));
    sigPosNeg = false(size(sigPosNeg));
end

return;




function iccp_ccpairs_crosscorr_strf_si_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)


% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);


% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);

hwPos = [ccpairs.hw_pos];

cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;



% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);

hwNeg = [ccpairs.hw_neg];

cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;


si = [ccpairs.strfsimilarity];


fprintf('\n');
fprintf('STRF Similarity population\n');
fprintf('N = %.0f\n', length(si));
fprintf('MN = %.4f\n', mean(si(~isnan(si))));
fprintf('SD = %.4f\n', std(si(~isnan(si) ) ) );
fprintf('MD = %.4f\n', median(si(~isnan(si))));
fprintf('MAD = %.4f\n', madstat(si(~isnan(si))));
fprintf('MX = %.4f\n', max(si(~isnan(si))) );


fprintf('\n');
fprintf('STRF Similarity: Pos Only\n');
siPos = si(sigPosOnly);
fprintf('N = %.0f\n', length(siPos));
fprintf('MN = %.4f\n', mean(siPos));
fprintf('SD = %.4f\n', std(siPos));
fprintf('MD = %.4f\n', median(siPos));
fprintf('MAD = %.4f\n', madstat(siPos));
fprintf('MX = %.4f\n', max(siPos));


fprintf('\n');
fprintf('STRF Similarity: Pos Neg Only\n');
siPosNeg = si(sigPosNeg);
fprintf('N = %.0f\n', length(siPosNeg));
fprintf('MN = %.4f\n', mean(siPosNeg));
fprintf('SD = %.4f\n', std(siPosNeg));
fprintf('MD = %.4f\n', median(siPosNeg));
fprintf('MAD = %.4f\n', madstat(siPosNeg));
fprintf('MX = %.4f\n', max(siPosNeg));


fprintf('\n');
fprintf('STRF Similarity: Neg Only\n');
siNeg = si(sigNegOnly);
fprintf('N = %.0f\n', length(siNeg));
fprintf('MN = %.4f\n', mean(siNeg));
fprintf('SD = %.4f\n', std(siNeg));
fprintf('MD = %.4f\n', median(siNeg));
fprintf('MAD = %.4f\n', madstat(siNeg));
fprintf('MX = %.4f\n', max(siNeg));



pPos_vs_PosNeg = ranksum(siPos, siPosNeg);
pPos_vs_Neg = ranksum(siPos, siNeg);
pPosNeg_vs_Neg = ranksum(siPosNeg, siNeg);

fprintf('\n');
fprintf('SI: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('SI: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('SI: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


fprintf('\n');

[r,p] = corrcoef(si(sigPosOnly), log10(cccPos(sigPosOnly)) );
fprintf('\nCCPairs: log10 Pos Only CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigPosOnly) )), log10(abs(cccPos(sigPosOnly) )));
fprintf('CCPairs: log10 Pos Only CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigPosOnly), log10(cccPos(sigPosOnly) ));
fprintf('CCPairs: log10 Pos Only CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

fprintf('\n');

[r,p] = corrcoef(si(sigPosNeg), log10(cccNeg(sigPosNeg)) );
fprintf('\nCCPairs: log10 Pos Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(log10(abs( si(sigPosNeg) )), log10(abs(cccNeg(sigPosNeg) )));
fprintf('CCPairs: log10 Pos Neg CCC vs. log10 SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));

[r,p] = corrcoef(si(sigPosNeg), log10(cccNeg(sigPosNeg) ));
fprintf('CCPairs: log10 Pos Neg CCC vs. SI: r = %.3f, p = %.3f\n', r(1,2), p(1,2));


fprintf('\n\n');


return;




function iccp_ccpairs_crosscorr_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)


% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);


% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);

hwPos = [ccpairs.hw_pos];

cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;



% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);

hwNeg = [ccpairs.hw_neg];

cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;



% index_pos_only = sigPos & ~sigNeg;
% index_neg_only = ~sigPos & sigNeg;
% index_pos_neg = sigPos & sigNeg;
% index_any = sigPos | sigNeg;
% 
% fprintf('Positive only peaks: %.0f\n', sum(index_pos_only) );
% fprintf('Negative only peaks: %.0f\n', sum(index_neg_only) );
% fprintf('Positive and Negative peaks: %.0f\n', sum(index_pos_neg) );
% fprintf('Any peaks: %.0f\n', sum(index_any) );


f = simple_stats(pdPos(sigPosOnly)); % facilitative only

s = simple_stats(pdNeg(sigPosNeg)); % suppressive

label = 'Facilitative';
fprintf('%s Peak Delay:\n', label);
fprintf('mean = %.3f\n', f.mn);
fprintf('SD = %.3f\n', f.sd);
fprintf('SE = %.3f\n', f.se);
fprintf('median = %.3f\n', f.md);
fprintf('MAD = %.3f\n', f.mad);

fprintf('\n');
label = 'Suppressive';
fprintf('%s Peak Delay:\n', label);
fprintf('mean = %.3f\n', s.mn);
fprintf('SD = %.3f\n', s.sd);
fprintf('SE = %.3f\n', s.se);
fprintf('median = %.3f\n', s.md);
fprintf('MAD = %.3f\n', s.mad);



f = simple_stats(hwPos(sigPosOnly));
s = simple_stats(hwNeg(sigPosNeg));

fprintf('\n');
label = 'Facilitative';
fprintf('%s Half-Width:\n', label);
fprintf('mean = %.3f\n', f.mn);
fprintf('SD = %.3f\n', f.sd);
fprintf('SE = %.3f\n', f.se);
fprintf('median = %.3f\n', f.md);
fprintf('MAD = %.3f\n', f.mad);
fprintf('\n');
label = 'Suppressive';
fprintf('%s Half-Width:\n', label);
fprintf('mean = %.3f\n', s.mn);
fprintf('SD = %.3f\n', s.sd);
fprintf('SE = %.3f\n', s.se);
fprintf('median = %.3f\n', s.md);
fprintf('MAD = %.3f\n', s.mad);


f = simple_stats(cccPos(sigPosOnly));
s = simple_stats(cccNeg(sigNegOnly));

fprintf('\n');
label = 'Facilitative';
fprintf('%s CCC:\n', label);
fprintf('mean = %.3f\n', f.mn);
fprintf('SD = %.3f\n', f.sd);
fprintf('SE = %.3f\n', f.se);
fprintf('median = %.3f\n', f.md);
fprintf('MAD = %.3f\n', f.mad);
fprintf('\n');
label = 'Suppressive';
fprintf('%s CCC:\n', label);
fprintf('mean = %.3f\n', s.mn);
fprintf('SD = %.3f\n', s.sd);
fprintf('SE = %.3f\n', s.se);
fprintf('median = %.3f\n', s.md);
fprintf('MAD = %.3f\n', s.mad);


return;




function iccp_sigfeature_to_fr_pli_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)

% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);

% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);
hwPos = [ccpairs.hw_pos];
cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;

% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);
hwNeg = [ccpairs.hw_neg];
cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;


fr1 = [ccpairs.fr1];
fr2 = [ccpairs.fr2];
fr = [fr1(:) fr2(:)];

pli1 = [ccpairs.pli1];
pli2 = [ccpairs.pli2];
pli = [pli1(:) pli2(:)];


frPos = fr(sigPosOnly, :);
frPosNeg = fr(sigPosNeg, :);
frNeg = fr(sigNegOnly, :);

pliPos = pli(sigPosOnly, :);
pliPosNeg = pli(sigPosNeg, :);
pliNeg = pli(sigNegOnly, :);



[larger,smaller] = iccp_largersmaller(fr(:,1), fr(:,2));
frdiffAll = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('\n');
fprintf('ALL FR vs FR: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');


[larger,smaller] = iccp_largersmaller(frPos(:,1), frPos(:,2));
frdiffPos = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('\n');
fprintf('Pos Only FR vs FR: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');


[larger,smaller] = iccp_largersmaller(frPosNeg(:,1), frPosNeg(:,2));
frdiffPosNeg = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('\n');
fprintf('Pos Neg FR vs FR: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');


[larger,smaller] = iccp_largersmaller(frNeg(:,1), frNeg(:,2));
frdiffNeg = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('\n');
fprintf('Neg Only FR vs FR: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');



fprintf('\n');
fprintf('All FR Difference\n');
fprintf('N = %.0f\n', length(frdiffAll));
fprintf('MN = %.4f\n', mean(frdiffAll));
fprintf('SD = %.4f\n', std(frdiffAll));
fprintf('MD = %.4f\n', median(frdiffAll));
fprintf('MAD = %.4f\n', madstat(frdiffAll));
fprintf('MX = %.4f\n', max(frdiffAll));


fprintf('\n');
fprintf('Pos Only FR Difference\n');
fprintf('N = %.0f\n', length(frdiffPos));
fprintf('MN = %.4f\n', mean(frdiffPos));
fprintf('SD = %.4f\n', std(frdiffPos));
fprintf('MD = %.4f\n', median(frdiffPos));
fprintf('MAD = %.4f\n', madstat(frdiffPos));
fprintf('MX = %.4f\n', max(frdiffPos));


fprintf('\n');
fprintf('Pos+Neg FR Difference\n');
fprintf('N = %.0f\n', length(frdiffPosNeg));
fprintf('MN = %.4f\n', mean(frdiffPosNeg));
fprintf('SD = %.4f\n', std(frdiffPosNeg));
fprintf('MD = %.4f\n', median(frdiffPosNeg));
fprintf('MAD = %.4f\n', madstat(frdiffPosNeg));
fprintf('MX = %.4f\n', max(frdiffPosNeg));


fprintf('\n');
fprintf('Neg Only FR Difference\n');
fprintf('N = %.0f\n', length(frdiffNeg));
fprintf('MN = %.4f\n', mean(frdiffNeg));
fprintf('SD = %.4f\n', std(frdiffNeg));
fprintf('MD = %.4f\n', median(frdiffNeg));
fprintf('MAD = %.4f\n', madstat(frdiffNeg));
fprintf('MX = %.4f\n', max(frdiffNeg));


pPos_vs_PosNeg = ranksum(frdiffPos, frdiffPosNeg);
pPos_vs_Neg = ranksum(frdiffPos, frdiffNeg);
pPosNeg_vs_Neg = ranksum(frdiffPosNeg, frdiffNeg);

fprintf('\n');
fprintf('FR diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('FR diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('FR diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');



[r,p] = corrcoef( log10(frdiffAll(sigPosOnly)), log10(cccPos(sigPosOnly)) );
fprintf('\n');
fprintf('EP: CCC vs FR Diff: r=%.4f, p = %.4f\n', r(2),p(2));
fprintf('\n');

[r,p] = corrcoef( log10(frdiffAll(sigPosNeg)), log10(cccNeg(sigPosNeg)) );
fprintf('\n');
fprintf('SP: CCC vs FR Diff: r=%.4f, p = %.4f\n', r(2),p(2));
fprintf('\n');




fprintf('\n\n');

[larger,smaller] = iccp_largersmaller(pli(:,1), pli(:,2));
plidiffAll = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL PLI vs PLI: r=%.4f, p = %.4f\n', r(2), p(2));


[larger,smaller] = iccp_largersmaller(pliPos(:,1), pliPos(:,2));
plidiffPos = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('Pos Only PLI vs PLI: r=%.4f, p = %.4f\n', r(2), p(2));


[larger,smaller] = iccp_largersmaller(pliPosNeg(:,1), pliPosNeg(:,2));
plidiffPosNeg = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('Pos Neg PLI vs PLI: r=%.4f, p = %.4f\n', r(2), p(2));


[larger,smaller] = iccp_largersmaller(pliNeg(:,1), pliNeg(:,2));
plidiffNeg = abs( larger - smaller );
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('Neg Only PLI vs PLI: r=%.4f, p = %.4f\n', r(2), p(2));



fprintf('\n');
fprintf('All PLI Difference\n');
fprintf('N = %.0f\n', length(plidiffAll));
fprintf('MN = %.4f\n', mean(plidiffAll));
fprintf('SD = %.4f\n', std(plidiffAll));
fprintf('MD = %.4f\n', median(plidiffAll));
fprintf('MAD = %.4f\n', madstat(plidiffAll));
fprintf('MX = %.4f\n', max(plidiffAll));


fprintf('\n');
fprintf('Pos Only PLI Difference\n');
fprintf('N = %.0f\n', length(plidiffPos));
fprintf('MN = %.4f\n', mean(plidiffPos));
fprintf('SD = %.4f\n', std(plidiffPos));
fprintf('MD = %.4f\n', median(plidiffPos));
fprintf('MAD = %.4f\n', madstat(plidiffPos));
fprintf('MX = %.4f\n', max(plidiffPos));


fprintf('\n');
fprintf('Pos+Neg PLI Difference\n');
fprintf('N = %.0f\n', length(plidiffPosNeg));
fprintf('MN = %.4f\n', mean(plidiffPosNeg));
fprintf('SD = %.4f\n', std(plidiffPosNeg));
fprintf('MD = %.4f\n', median(plidiffPosNeg));
fprintf('MAD = %.4f\n', madstat(plidiffPosNeg));
fprintf('MX = %.4f\n', max(plidiffPosNeg));


fprintf('\n');
fprintf('Neg Only PLI Difference\n');
fprintf('N = %.0f\n', length(plidiffNeg));
fprintf('MN = %.4f\n', mean(plidiffNeg));
fprintf('SD = %.4f\n', std(plidiffNeg));
fprintf('MD = %.4f\n', median(plidiffNeg));
fprintf('MAD = %.4f\n', madstat(plidiffNeg));
fprintf('MX = %.4f\n', max(plidiffNeg));


pPos_vs_PosNeg = ranksum(plidiffPos, plidiffPosNeg);
pPos_vs_Neg = ranksum(plidiffPos, plidiffNeg);
pPosNeg_vs_Neg = ranksum(plidiffPosNeg, plidiffNeg);

fprintf('\n');
fprintf('PLI diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('PLI diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('PLI diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');



[r,p] = corrcoef( log10(plidiffAll(sigPosOnly)), log10(cccPos(sigPosOnly)) );
fprintf('\n');
fprintf('EP: CCC vs FR Diff: r=%.4f, p = %.4f\n', r(2),p(2));
fprintf('\n');

[r,p] = corrcoef( log10(plidiffAll(sigPosNeg)), log10(cccNeg(sigPosNeg)) );
fprintf('\n');
fprintf('SP: CCC vs FR Diff: r=%.4f, p = %.4f\n', r(2),p(2));
fprintf('\n');

return;




function iccp_sigfeature_to_cf_q_lat_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)

% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);

% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);
hwPos = [ccpairs.hw_pos];
cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;

% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);
hwNeg = [ccpairs.hw_neg];
cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;

% CF
cf1 = [ccpairs.cf1];
cf2 = [ccpairs.cf2];
cf = [cf1(:) cf2(:)];
cfPos = cf(sigPosOnly,:);
cfPosNeg = cf(sigPosNeg,:);
cfNeg = cf(sigNegOnly,:);

% Q - spectral tuning
q1 = [ccpairs.q1];
q2 = [ccpairs.q2];
q = [q1(:) q2(:)];
qPos = q(sigPosOnly,:);
qPosNeg = q(sigPosNeg,:);
qNeg = q(sigNegOnly,:);

% Latency
latency1 = [ccpairs.latency1];
latency2 = [ccpairs.latency2];
latency = [latency1(:) latency2(:)];
latencyPos = latency(sigPosOnly,:);
latencyPosNeg = latency(sigPosNeg,:);
latencyNeg = latency(sigNegOnly,:);


[larger,smaller] = iccp_largersmaller(cf(:,1), cf(:,2));
cfDiffAll = abs( log2(larger./smaller) );
cfDiffPosOnly = cfDiffAll(sigPosOnly);
cfDiffPosNeg = cfDiffAll(sigPosNeg);
cfDiffNegOnly = cfDiffAll(sigNegOnly);

fprintf('\n');
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL CF vs CF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosOnly)), log10(smaller(sigPosOnly)) );
fprintf('Pos Only CF vs CF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosNeg)), log10(smaller(sigPosNeg)) );
fprintf('Pos+Neg CF vs CF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigNegOnly)), log10(smaller(sigNegOnly)) );
fprintf('Neg Only CF vs CF: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');



[larger,smaller] = iccp_largersmaller(q(:,1), q(:,2));
qDiffAll = abs( log2(larger./smaller) );
qDiffPosOnly = qDiffAll(sigPosOnly);
qDiffPosNeg = qDiffAll(sigPosNeg);
qDiffNegOnly = qDiffAll(sigNegOnly);

fprintf('\n');
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL Q vs Q: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosOnly)), log10(smaller(sigPosOnly)) );
fprintf('Pos Only Q vs Q: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosNeg)), log10(smaller(sigPosNeg)) );
fprintf('Pos+Neg Q vs Q: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigNegOnly)), log10(smaller(sigNegOnly)) );
fprintf('Neg Only Q vs Q: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');




[larger,smaller] = iccp_largersmaller(latency(:,1), latency(:,2));
latDiffAll = abs( larger - smaller );
latDiffPosOnly = latDiffAll(sigPosOnly);
latDiffPosNeg = latDiffAll(sigPosNeg);
latDiffNegOnly = latDiffAll(sigNegOnly);

fprintf('\n');
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL Lat vs Lat: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosOnly)), log10(smaller(sigPosOnly)) );
fprintf('Pos Only Lat vs Lat: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosNeg)), log10(smaller(sigPosNeg)) );
fprintf('Pos+Neg Lat vs Lat: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigNegOnly)), log10(smaller(sigNegOnly)) );
fprintf('Neg Only Lat vs Lat: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');



% CF difference analysis

fprintf('\n');
fprintf('All CF Difference\n');
data = cfDiffAll;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Pos Only CF Difference\n');
data = cfDiffPosOnly;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg CF Difference\n');
data = cfDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only CF Difference\n');
data = cfDiffNegOnly;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(cfDiffPosOnly, cfDiffPosNeg);
pPos_vs_Neg = ranksum(cfDiffPosOnly, cfDiffNegOnly);
pPosNeg_vs_Neg = ranksum(cfDiffPosNeg, cfDiffNegOnly);

fprintf('\n');
fprintf('CF diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('CF diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('CF diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');

[r,p] = corrcoef( log10(cfDiffAll(sigPosOnly)), log10(cccPos(sigPosOnly)) );
fprintf('EP: CCC vs CF Diff: r=%.4f, p = %.4f\n', r(2),p(2));

[r,p] = corrcoef( log10(cfDiffAll(sigPosNeg)), log10(cccNeg(sigPosNeg)) );
fprintf('SP: CCC vs CF Diff: r=%.4f, p = %.4f\n', r(2),p(2));




% Q difference analysis

fprintf('\n');
fprintf('All Q Difference\n');
data = qDiffAll;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos Only Q Difference\n');
data = qDiffPosOnly;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg Q Difference\n');
data = qDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only Q Difference\n');
data = qDiffNegOnly;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(qDiffPosOnly, qDiffPosNeg);
pPos_vs_Neg = ranksum(qDiffPosOnly, qDiffNegOnly);
pPosNeg_vs_Neg = ranksum(qDiffPosNeg, qDiffNegOnly);

fprintf('\n');
fprintf('Q diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('Q diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('Q diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');

[r,p] = corrcoef( log10(qDiffAll(sigPosOnly)), log10(cccPos(sigPosOnly)) );
fprintf('EP: CCC vs Q Diff: r=%.4f, p = %.4f\n', r(2),p(2));

[r,p] = corrcoef( log10(qDiffAll(sigPosNeg)), log10(cccNeg(sigPosNeg)) );
fprintf('SP: CCC vs Q Diff: r=%.4f, p = %.4f\n', r(2),p(2));



% Latency difference analysis

fprintf('\n');
fprintf('All Latency Difference\n');
data = latDiffAll;
simple_stats(data,1);


fprintf('\n');
fprintf('Pos Only Latency Difference\n');
data = latDiffPosOnly;
simple_stats(data,1);


fprintf('\n');
fprintf('Pos+Neg Latency Difference\n');
data = latDiffPosNeg;
simple_stats(data,1);


fprintf('\n');
fprintf('Neg Only Latency Difference\n');
data = latDiffNegOnly;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(latDiffPosOnly, latDiffPosNeg);
pPos_vs_Neg = ranksum(latDiffPosOnly, latDiffNegOnly);
pPosNeg_vs_Neg = ranksum(latDiffPosNeg, latDiffNegOnly);

fprintf('\n');
fprintf('Latency diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('Latency diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('Latency diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');

[r,p] = corrcoef( log10(latDiffAll(sigPosOnly)), log10(cccPos(sigPosOnly)) );
fprintf('EP: CCC vs Latency Diff: r=%.4f, p = %.4f\n', r(2),p(2));

[r,p] = corrcoef( log10(latDiffAll(sigPosNeg)), log10(cccNeg(sigPosNeg)) );
fprintf('SP: CCC vs Latency Diff: r=%.4f, p = %.4f\n', r(2),p(2));

return;




function iccp_sigfeature_to_bmf_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)

% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);

% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);
hwPos = [ccpairs.hw_pos];
cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;

% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);
hwNeg = [ccpairs.hw_neg];
cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;

% bTMF
btmf1 = [ccpairs.btmf1];
btmf2 = [ccpairs.btmf2];
btmf = [btmf1(:) btmf2(:)];

btmfPos = btmf(sigPosOnly,:);
btmfPosNeg = btmf(sigPosNeg,:);
btmfNeg = btmf(sigNegOnly,:);


% bSMF
bsmf1 = [ccpairs.bsmf1];
bsmf2 = [ccpairs.bsmf2];
bsmf = [bsmf1(:) bsmf2(:)];

bsmfPos = bsmf(sigPosOnly,:);
bsmfPosNeg = bsmf(sigPosNeg,:);
bsmfNeg = bsmf(sigNegOnly,:);




% bTMF differences
[larger,smaller] = iccp_largersmaller(btmf(:,1), btmf(:,2));
btmfDiffAll = abs( larger - smaller );
btmfDiffPosOnly = btmfDiffAll(sigPosOnly);
btmfDiffPosNeg = btmfDiffAll(sigPosNeg);
btmfDiffNegOnly = btmfDiffAll(sigNegOnly);

fprintf('\n');
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL bTMF vs bTMF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosOnly)), log10(smaller(sigPosOnly)) );
fprintf('Pos Only bTMF vs bTMF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosNeg)), log10(smaller(sigPosNeg)) );
fprintf('Pos+Neg bTMF vs bTMF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigNegOnly)), log10(smaller(sigNegOnly)) );
fprintf('Neg Only bTMF vs bTMF: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');


% bTMF difference statistics
fprintf('\n');
fprintf('All BTMF Difference\n');
data = btmfDiffAll;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Pos Only bTMF Difference\n');
data = btmfDiffPosOnly;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg bTMF Difference\n');
data = btmfDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only bTMF Difference\n');
data = btmfDiffNegOnly;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(btmfDiffPosOnly, btmfDiffPosNeg);
pPos_vs_Neg = ranksum(btmfDiffPosOnly, btmfDiffNegOnly);
pPosNeg_vs_Neg = ranksum(btmfDiffPosNeg, btmfDiffNegOnly);

fprintf('\n');
fprintf('bTMF diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('bTMF diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('bTMF diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');

x = log10(btmfDiffAll(sigPosOnly));
y = log10(cccPos(sigPosOnly));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('EP: CCC vs bTMF Diff: r=%.4f, p = %.4f\n', r(2),p(2));


x = log10(btmfDiffAll(sigPosNeg));
y = log10(cccNeg(sigPosNeg));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('SP: CCC vs bTMF Diff: r=%.4f, p = %.4f\n', r(2),p(2));




% bSMF differences
[larger,smaller] = iccp_largersmaller(bsmf(:,1), bsmf(:,2));
bsmfDiffAll = abs( larger - smaller );
bsmfDiffPosOnly = bsmfDiffAll(sigPosOnly);
bsmfDiffPosNeg = bsmfDiffAll(sigPosNeg);
bsmfDiffNegOnly = bsmfDiffAll(sigNegOnly);

fprintf('\n');
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL bSMF vs bSMF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosOnly)), log10(smaller(sigPosOnly)) );
fprintf('Pos Only bSMF vs bSMF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosNeg)), log10(smaller(sigPosNeg)) );
fprintf('Pos+Neg bSMF vs bSMF: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigNegOnly)), log10(smaller(sigNegOnly)) );
fprintf('Neg Only bSMF vs bSMF: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');


% bSMF difference statistics
fprintf('\n');
fprintf('All bSMF Difference\n');
data = bsmfDiffAll;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Pos Only bSMF Difference\n');
data = bsmfDiffPosOnly;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg bSMF Difference\n');
data = bsmfDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only bSMF Difference\n');
data = bsmfDiffNegOnly;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(bsmfDiffPosOnly, bsmfDiffPosNeg);
pPos_vs_Neg = ranksum(bsmfDiffPosOnly, bsmfDiffNegOnly);
pPosNeg_vs_Neg = ranksum(bsmfDiffPosNeg, bsmfDiffNegOnly);

fprintf('\n');
fprintf('bSMF diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('bSMF diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('bSMF diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


x = log10(bsmfDiffAll(sigPosOnly));
y = log10(cccPos(sigPosOnly));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('EP: CCC vs bSMF Diff: r=%.4f, p = %.4f\n', r(2),p(2));


x = log10(bsmfDiffAll(sigPosNeg));
y = log10(cccNeg(sigPosNeg));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('SP: CCC vs bSMF Diff: r=%.4f, p = %.4f\n', r(2),p(2));

return;




function iccp_sigfeature_to_fio_si_asi_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)

% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);

% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);
hwPos = [ccpairs.hw_pos];
cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;

% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);
hwNeg = [ccpairs.hw_neg];
cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;


% Get ASI nonlinearity metric
asi1 = [ccpairs.fioasi1];
asi2 = [ccpairs.fioasi2];
asi = [asi1(:) asi2(:)];
asiPos = asi(sigPosOnly,:);
asiPosNeg = asi(sigPosNeg,:);
asiNeg = asi(sigNegOnly,:);


% Get SI nonlinearity metric
si = [ccpairs.fiosimilarity];
siPos = si(sigPosOnly);
siPosNeg = si(sigPosNeg);
siNeg = si(sigNegOnly);




% ASI differences
[larger,smaller] = iccp_largersmaller(asi(:,1), asi(:,2));
asiDiffAll = abs( larger - smaller );
asiDiffPosOnly = asiDiffAll(sigPosOnly);
asiDiffPosNeg = asiDiffAll(sigPosNeg);
asiDiffNegOnly = asiDiffAll(sigNegOnly);

fprintf('\n');
[r,p] = corrcoef( log10(larger), log10(smaller) );
fprintf('ALL ASI vs ASI: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosOnly)), log10(smaller(sigPosOnly)) );
fprintf('Pos Only ASI vs ASI: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigPosNeg)), log10(smaller(sigPosNeg)) );
fprintf('Pos+Neg ASI vs ASI: r=%.4f, p = %.4f\n', r(2), p(2));

[r,p] = corrcoef( log10(larger(sigNegOnly)), log10(smaller(sigNegOnly)) );
fprintf('Neg Only ASI vs ASI: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');


% ASI difference statistics
fprintf('\n');
fprintf('All ASI Difference\n');
data = asiDiffAll;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Pos Only ASI Difference\n');
data = asiDiffPosOnly;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg ASI Difference\n');
data = asiDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only ASI Difference\n');
data = asiDiffNegOnly;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(asiDiffPosOnly, asiDiffPosNeg);
pPos_vs_Neg = ranksum(asiDiffPosOnly, asiDiffNegOnly);
pPosNeg_vs_Neg = ranksum(asiDiffPosNeg, asiDiffNegOnly);

fprintf('\n');
fprintf('ASI diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('ASI diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('ASI diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


x = log10(asiDiffAll(sigPosOnly));
y = log10(cccPos(sigPosOnly));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('EP: CCC vs ASI Diff: r=%.4f, p = %.4f\n', r(2),p(2));


x = log10(asiDiffAll(sigPosNeg));
y = log10(cccNeg(sigPosNeg));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('SP: CCC vs ASI Diff: r=%.4f, p = %.4f\n', r(2),p(2));




% SI difference statistics
fprintf('\n');
fprintf('All FIO SI\n');
data = si;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Pos Only FIO SI\n');
data = siPos;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Pos+Neg FIO SI\n');
data = siPosNeg;
simple_stats(data,1);
fprintf('\n');

fprintf('\n');
fprintf('Neg Only FIO SI\n');
data = siNeg;
simple_stats(data,1);
fprintf('\n');


pPos_vs_PosNeg = ranksum(siPos, siPosNeg);
pPos_vs_Neg = ranksum(siPos, siNeg);
pPosNeg_vs_Neg = ranksum(siPosNeg, siNeg);

fprintf('\n');
fprintf('FIO SI: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('FIO SI: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('FIO SI: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


x = si(sigPosOnly);
y = log10(cccPos(sigPosOnly));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('EP: CCC vs SI: r=%.4f, p = %.4f\n', r(2),p(2));


x = si(sigPosNeg);
y = log10(cccNeg(sigPosNeg));
index = ~isinf(x(:)) & ~isinf(y(:));
[r,p] = corrcoef( x(index), y(index) );
fprintf('SP: CCC vs SI: r=%.4f, p = %.4f\n', r(2),p(2));



return;




function iccp_sigfeature_to_fiofit_stats(ccpairs, sigPosOnly, sigPosNeg, sigNegOnly)

% Any cross-cov function that had a positive peak
sigPos = [ccpairs.sigfeature_pos];
sigPos = logical(sigPos);

% Any cross-cov function that had a positive peak
sigNeg = [ccpairs.sigfeature_neg];
sigNeg = logical(sigNeg);

% Positive peaks
pdPos = [ccpairs.pd_pos];
pdPos = abs(pdPos);
hwPos = [ccpairs.hw_pos];
cccPos = [ccpairs.ccc_pos];
cccPos(cccPos < 0) = 0.0001;

% Negative peaks
pdNeg = [ccpairs.pd_neg];
pdNeg = abs(pdNeg);
hwNeg = [ccpairs.hw_neg];
cccNeg = [ccpairs.ccc_neg];
cccNeg(cccNeg < 0) = 0.0001;


% Get ASI nonlinearity metric
asi1 = [ccpairs.fioasi1];
asi2 = [ccpairs.fioasi2];
asi = [asi1(:) asi2(:)];
asiPos = asi(sigPosOnly,:);
asiPosNeg = asi(sigPosNeg,:);
asiNeg = asi(sigNegOnly,:);


% Get SI nonlinearity metric
si = [ccpairs.fiosimilarity];
siPos = si(sigPosOnly);
siPosNeg = si(sigPosNeg);
siNeg = si(sigNegOnly);



% Get nonlinearity metrics
theta1 = [ccpairs.theta1];
theta2 = [ccpairs.theta2];
theta = [theta1(:) theta2(:)];



sigma1 = [ccpairs.sigma1];
sigma2 = [ccpairs.sigma2];
sigma = [sigma1(:) sigma2(:)];



nmse1 = [ccpairs.nmse1];
index_sig = find(nmse1 <= 0.2);
index_nonsig = find(nmse1 > 0.2);
nmse1(index_sig) = 1;
nmse1(index_nonsig) = 0;
nmse1 = logical(nmse1);

nmse2 = [ccpairs.nmse2];
index_sig = find(nmse2 <= 0.2);
index_nonsig = find(nmse2 > 0.2);
nmse2(index_sig) = 1;
nmse2(index_nonsig) = 0;
nmse2 = logical(nmse2);

nmse = nmse1 & nmse2;


r21 = [ccpairs.r21];
index_sig = find(r21 >=0.8);
index_nonsig = find(r21 < 0.8);
r21(index_sig) = 1;
r21(index_nonsig) = 0;
r21 = logical(r21);

r22 = [ccpairs.r22];
index_sig = find(r22 >=0.8);
index_nonsig = find(r22 < 0.8);
r22(index_sig) = 1;
r22(index_nonsig) = 0;
r22 = logical(r22);

r2 = r21 & r22;


% Theta analysis
[larger, smaller] = iccp_largersmaller(theta(:,1), theta(:,2));


fprintf('\n');
indexAll = nmse & r2;
[r,p] = corrcoef( larger(indexAll), smaller(indexAll) );
fprintf('ALL FIOFIT Theta vs Theta: r=%.4f, p = %.4f\n', r(2), p(2));


indexPos = nmse & r2 & sigPosOnly;
[r,p] = corrcoef( larger(indexPos), smaller(indexPos) );
fprintf('Pos Only FIOFIT Theta vs Theta: r=%.4f, p = %.4f\n', r(2), p(2));


indexPosNeg = nmse & r2 & sigPosNeg;
[r,p] = corrcoef( larger(indexPosNeg), smaller(indexPosNeg) );
fprintf('Pos+Neg FIOFIT Theta vs Theta: r=%.4f, p = %.4f\n', r(2), p(2));


indexNeg = nmse & r2 & sigNegOnly;
[r,p] = corrcoef( larger(indexNeg), smaller(indexNeg) );
fprintf('Neg Only FIOFIT Theta vs Theta: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');



% Difference distributions
thetaDiffAll = abs( larger(indexAll) - smaller(indexAll) );
thetaDiffPos = abs( larger(indexPos) - smaller(indexPos) );
thetaDiffPosNeg = abs( larger(indexPosNeg) - smaller(indexPosNeg) );
thetaDiffNeg = abs( larger(indexNeg) - smaller(indexNeg) );



% Difference statistics
fprintf('\n');
fprintf('All theta Difference\n');
data = thetaDiffAll;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos Only theta Difference\n');
data = thetaDiffPos;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg theta Difference\n');
data = thetaDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only theta Difference\n');
data = thetaDiffNeg;
simple_stats(data,1);



pPos_vs_PosNeg = ranksum(thetaDiffPos, thetaDiffPosNeg);
pPos_vs_Neg = ranksum(thetaDiffPos, thetaDiffNeg);
pPosNeg_vs_Neg = ranksum(thetaDiffPosNeg, thetaDiffNeg);

fprintf('\n');
fprintf('FIOFIT Theta Diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('FIOFIT Theta Diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('FIOFIT Theta Diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


% x = si(sigPosOnly);
% y = log10(cccPos(sigPosOnly));
% index = ~isinf(x(:)) & ~isinf(y(:));
% [r,p] = corrcoef( x(index), y(index) );
% fprintf('EP: CCC vs SI: r=%.4f, p = %.4f\n', r(2),p(2));
% 
% 
% x = si(sigPosNeg);
% y = log10(cccNeg(sigPosNeg));
% index = ~isinf(x(:)) & ~isinf(y(:));
% [r,p] = corrcoef( x(index), y(index) );
% fprintf('SP: CCC vs SI: r=%.4f, p = %.4f\n', r(2),p(2));




% Sigma analysis
[larger, smaller] = iccp_largersmaller(sigma(:,1), sigma(:,2));


fprintf('\n');
indexAll = nmse & r2;
[r,p] = corrcoef( larger(indexAll), smaller(indexAll) );
fprintf('ALL FIOFIT Sigma vs Sigma: r=%.4f, p = %.4f\n', r(2), p(2));


indexPos = nmse & r2 & sigPosOnly;
[r,p] = corrcoef( larger(indexPos), smaller(indexPos) );
fprintf('Pos Only FIOFIT Sigma vs Sigma: r=%.4f, p = %.4f\n', r(2), p(2));


indexPosNeg = nmse & r2 & sigPosNeg;
[r,p] = corrcoef( larger(indexPosNeg), smaller(indexPosNeg) );
fprintf('Pos+Neg FIOFIT Sigma vs Sigma: r=%.4f, p = %.4f\n', r(2), p(2));


indexNeg = nmse & r2 & sigNegOnly;
[r,p] = corrcoef( larger(indexNeg), smaller(indexNeg) );
fprintf('Neg Only FIOFIT Sigma vs Sigma: r=%.4f, p = %.4f\n', r(2), p(2));
fprintf('\n');



% Difference distributions
sigmaDiffAll = abs( larger(indexAll) - smaller(indexAll) );
sigmaDiffPos = abs( larger(indexPos) - smaller(indexPos) );
sigmaDiffPosNeg = abs( larger(indexPosNeg) - smaller(indexPosNeg) );
sigmaDiffNeg = abs( larger(indexNeg) - smaller(indexNeg) );



% Difference statistics
fprintf('\n');
fprintf('All sigma Difference\n');
data = sigmaDiffAll;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos Only sigma Difference\n');
data = sigmaDiffPos;
simple_stats(data,1);

fprintf('\n');
fprintf('Pos+Neg sigma Difference\n');
data = sigmaDiffPosNeg;
simple_stats(data,1);

fprintf('\n');
fprintf('Neg Only sigma Difference\n');
data = sigmaDiffNeg;
simple_stats(data,1);


pPos_vs_PosNeg = ranksum(sigmaDiffPos, sigmaDiffPosNeg);
pPos_vs_Neg = ranksum(sigmaDiffPos, sigmaDiffNeg);
pPosNeg_vs_Neg = ranksum(sigmaDiffPosNeg, sigmaDiffNeg);

fprintf('\n');
fprintf('FIOFIT Sigma Diff: Pos vs PosNeg: p = %.4f\n', pPos_vs_PosNeg);
fprintf('FIOFIT Sigma Diff: Pos vs Neg: p = %.4f\n', pPos_vs_Neg);
fprintf('FIOFIT Sigma Diff: PosNeg vs Neg: p = %.4f\n', pPosNeg_vs_Neg);
fprintf('\n');


return;




