function index = iccpairs_find_fra_params_element(params, exp, site, chan, model, stim)
% iccpairs_find_fra_params_element Index to struct array element for experiment specs
% 
%     index = iccpairs_find_fra_params_element(params, exp, site, chan, model, stim)
% 
%     params : struct array holding STRF data
% 
%     exp : experiment of data that is sought
%     site : site of data
%     chan : channel of data
%     model : model of data
% 
%     index : the index into params so that params(index) contains the input
%     argument information.
% 
%     If index = 0, then element was not found.
% 
%     index = iccpairs_find_fra_params_element(params, exp, site, chan, model)

model = sort(model);
index = 0;

for i = 1:length(params)

    if ( strcmp(params(i).exp, exp) )
        if ( params(i).site == site )
            if ( params(i).chan == chan ) 
                if ( strcmp(params(i).stim, stim) ) 
                    params_model = sort(params(i).model);
                    if ( params_model(1) == model(1) )
                        index = i;
                    end
                end
            end
        end
    end

end % (for i)

return;

