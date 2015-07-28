function [fiodata] = iccp_ccpairs_select_max_nonlinearity(fiodata)
% iccpairs_ccpairs_select_max_plot_strf_params STRF params for synchronous spikes
%
% [fiodata] = iccp_ccpairs_select_max_plot_strf_ptparams;
% iccp_ccpairs_select_max_plot_strf_params(fiodata)
% --------------------------------------------------------
%     fradata : struct holding the pairwise FRA parameters for synchronous
%     spikes.
% 
%     The call:
% 
%         [rfdata] = iccpairs_ccpairs_select_max_plot_strf_params;
% 
%     Reads through the files in the current folder:
% 
%     *-strfcmb-pairs-params.mat     
%     *-strfcmb-pairs-isi-fixed-ccpairs-select-max.mat   
% 
%     and calculates every pairwise combination of frequency response
%     area parameters for each electrode channel. The results are stored and
%     returned in the struct fradata.
% 
%     The pairs are then compared to the struct array ccpairs that is 
%     stored in the ccpairs-select-max.mat file.
% 
%     If the neuron pair for the FRA param data is an element of ccpairs,
%     then the data is kept and later plotted. Data that do not have a match
%     in ccpairs is discarded.
% 
%     This procedure ensures that only data from synchronous spiking neurons
%     is plotted. This also ensures that only resolved neurons have their data
%     plotted. Only synchronous spikes in the correlogram guarantees that
%     the spike sorting procedure worked appropriately. This check implicitly
%     controls for this.
% 
%     The call:
% 
%         [rfdata] = iccpairs_ccpairs_select_max_plot_strf_params;
% 
%     produces the struct rfdata.
% 
%     The call:
% 
%         iccpairs_ccpairs_select_max_plot_fra_params(rfdata);
% 
%     uses the already created fradata struct and saves time, since this
%     skips reading through all the files in the current folder.
% 
%     iccpairs_ccpairs_select_max_plot_fra_params(rfdata);
%     [rfdata] = iccpairs_ccpairs_select_max_plot_strf_params;
%


nargchk(0,1,nargin);

if ( nargin == 0 )
    [position, cc, theta, sigma, asi] = get_pairs_fioparams;
    fiodata.cc = cc;
    fiodata.position = position;
    fiodata.theta = theta;
    fiodata.sigma = sigma;
    fiodata.asi = asi;
end

return;



function [position, cc, theta, sigma, asi] = get_pairs_fioparams

position = [];
cc = [];
theta = [];
sigma = [];
asi = [];

d = dir('*-strfcmb-pairs-isi-fixed-ccpairs-select-max.mat');

fiopath = 'D:\ICC_Pairs_Data\sta-fio-fitparams'

for n = 1:length(d)

    ccfile = d(n).name;
    s = load(ccfile, 'ccpairs_max');
    ccpairs = s.ccpairs_max;

    if ( length(ccpairs) > 0 )

        index = findstr(ccfile, '-isi-fixed-ccpairs-select-max.mat');
        prefix = ccfile(1:index-1);

        paramsfile = fullfile(fiopath, sprintf('%s-sta-fio-fitparams.mat', prefix));

        s = load(paramsfile, 'fioFit');
        fioFit = s.fioFit;

        for ii = 1:length(ccpairs)

            % specs for the neuron pair
            exp = ccpairs(ii).exp;
            site = ccpairs(ii).site;
            chan = ccpairs(ii).chan;
            stim = ccpairs(ii).stim;

            model1 = ccpairs(ii).model1; % spike model for neuron 1
            model2 = ccpairs(ii).model2; % for neuron 2


            % Find where the parameter values are for each neuron
            index1 = iccpairs_find_fra_params_element(fioFit, exp, ...
                site, chan, model1, stim);

            index2 = iccpairs_find_fra_params_element(fioFit, exp, ...
                site, chan, model2, stim);

            % Make sure we found the data
            if ( index1 > 0 && index2 > 0 )

                if fioFit(index1).nmse_sta < 0.1 && fioFit(index2).nmse_sta < 0.1
                    cc = [cc; ccpairs(ii).ccc];
                    position = [position; fioFit(index1).position];

                    sigma1 = fioFit(index1).fitParams_sta(3);
                    theta1 = fioFit(index1).fitParams_sta(2);
                    sigma2 = fioFit(index2).fitParams_sta(3);
                    theta2 = fioFit(index2).fitParams_sta(2);
                    asi1 = fioFit(index1).asi_sta;
                    asi2 = fioFit(index2).asi_sta;

                    sigma = [sigma; sigma1 sigma2];
                    theta = [theta; theta1 theta2];
                    asi = [asi; asi1 asi2];
                end

            end

        end % (for ii)

    end % (if)

end % (for n)

return;











