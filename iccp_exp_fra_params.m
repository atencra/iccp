function [center, range, nreps] = iccp_exp_fra_params(exp, site, depth)
% iccp_exp_fra_params TMS tuning curve parameters used in cat ICC experiments
% 
%     [center, range, nreps] = iccp_exp_fra_params(exp, site, depth)
% 
%     Given the experiment date, site, and michigan probe insertion depth, 
%     this function returns the TMS system paramters that were used to 
%     present pure tones. 
% 
%     The TMS parameters are the center frequency of the tones presented,
%     the range of the frequencies, in octaves, and the the number
%     of tone repeats. For the TMS system, if a tone was repeated,
%     then the tone was played repeated in sequence. Thus, if tone A
%     was repeated five times, followed by tone B, then we would have
%     the stimulus sequence:
% 
%     A A A A A B B B B B ...
% 
%     center : center frequency, in kHz
% 
%     range : of frequencies, in octaves
% 
%     nreps : integer number of repeats, 1-N




if ( strcmp(exp, '2004-2-17') )

    if ( site == 1 && depth == 2395 )
        center = 1.6;
        range = 6;
        nreps = 5;

    elseif ( site == 4 && depth == 2405 )
        center = 1.6;
        range = 6;
        nreps = 5;

    elseif ( site == 11 && depth == 4046 )
        center = 2.5;
        range = 6;
        nreps = 5;


    elseif ( site == 12 && depth == 4270 )
        center = 5;
        range = 6;
        nreps = 5;

    elseif ( site == 13 && depth == 4057 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 14 && depth == 5000 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 14 && depth == 5000 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 19 && depth == 4767 )
        center = 2.5;
        range = 6;
        nreps = 5;

    elseif ( site == 19 && depth == 5005 )
        center = 2.5;
        range = 6;
        nreps = 3;
    else
        error('No FRA parameters for the exp/site/depth.');
    end

elseif ( strcmp(exp,'2003-11-24') )

    if ( ismember(site, 1:6) )
        center = 1.6;
        range = 6;
        nreps = 5;

    elseif ( ismember(site, 7:13) )
        center = 1;
        range = 6;
        nreps = 5;

    elseif ( site == 14 )
        if ( depth == 2420 )
            center = 1;
            range = 6;
            nreps = 5;
        elseif ( depth == 4620 )
        else
            error('Bad depth.');
        end
    else
        error('No FRA parameters for the exp/site/depth.');
    end

end

return;
