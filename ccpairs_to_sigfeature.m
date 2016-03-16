function [pdPos, hwPos, cccPos, pdNeg, hwNeg, cccNeg, sigPos, sigNeg] = ...
    ccpairs_to_sigfeature(ccpairs)

if ( isfield(ccpairs,'pd_pos') && isfield(ccpairs,'pd_neg') )

    pdPos = [ccpairs.pd_pos];
    pdPos = abs(pdPos); % only absolute value matters, not direction
    hwPos = [ccpairs.hw_pos];
    sigPos = [ccpairs.sigfeature_pos];
    sigPos = logical(sigPos);
    cccPos = [ccpairs.ccc_pos];
    cccPos(cccPos < 0) = 0.0001;

    pdNeg = [ccpairs.pd_neg];
    pdNeg = abs(pdNeg);
    hwNeg = [ccpairs.hw_neg];
    sigNeg = [ccpairs.sigfeature_neg];
    sigNeg = logical(sigNeg & pdNeg > 0); %***
    cccNeg = [ccpairs.ccc_neg];
    cccNeg(cccNeg < 0) = 0.0001;

    index_pos_only = sigPos & ~sigNeg;
    index_neg_only = ~sigPos & sigNeg;
    index_pos_neg = sigPos & sigNeg;
    index_any = sigPos | sigNeg;

    fprintf('Positive only peaks: %.0f\n', sum(index_pos_only) );
    fprintf('Negative only peaks: %.0f\n', sum(index_neg_only) );
    fprintf('Positive and Negative peaks: %.0f\n', sum(index_pos_neg) );
    fprintf('Any peaks: %.0f\n', sum(index_any) );

end



if ( ~isfield(ccpairs,'pd_pos') && isfield(ccpairs,'peakdelay') )

    pdPos = [ccpairs.peakdelay];
    pdPos = abs(pdPos);
    hwPos = [ccpairs.halfwidth];
    sigPos = [ccpairs.significant];
    sigPos = logical(sigPos);
    cccPos = [ccpairs.ccc];
    cccPos(cccPos < 0) = 0.0001;

    pdNeg = zeros(size(sigPos)); 
    hwNeg = zeros(size(sigPos));
    cccNeg = zeros(size(sigPos));
    sigNeg = false(size(sigPos));

end


return;