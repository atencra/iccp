function [nb_nb_ne, ne_nb_ne] = iccp_plot_pairs_fra_resptype(nb_nb_ne, ne_nb_ne)
% iccp_plot_pairs_fra_resptype Compare response types of ICC pairs
%
% [nb_nb_ne, ne_nb_ne] = iccp_plot_pairs_fra_resptype(nb_nb_ne, ne_nb_ne)
% ---------------------------------------------------------------------
% Reads through files in a directory to find struct arrays
% holding the response type data. The files have names in the following
% form:
%        *-fracmb-pairs-resptype.mat
%
% The function finds the neurons that were recorded from the 
% same channel, and then compares the response type parameters
% for these neurons.
%
% Function calls:
%
% iccp_plot_pairs_fra_resptype : searches through files and plots the data
% 
% [nb_nb_ne, ne_nb_ne] = iccp_plot_pairs_fra_resptype : returns the data used
% to make the plots.
% 
% iccp_plot_pairs_fra_resptype(nb_nb_ne, ne_nb_ne) : uses the previously 
% returned data to make the plots. This speeds up the process considerably,
% since there could be many files.


% Get data from files if no input arguments
if ( nargin == 0 )
   [position, nb_nb_ne, ne_nb_ne] = get_pairs_fra_resptype;
end

close all;

% Plot response type data - scatter and diff histograms
iccp_plot_pairs_fra_resptype_scatter_hist(nb_nb_ne, ne_nb_ne);


% Randomization test for response types
iccp_plot_pairs_fra_resptype_diff_rand(nb_nb_ne, ne_nb_ne);

return;



function [position, nb_nb_ne, ne_nb_ne] = get_pairs_fra_resptype

position = [];
nb_nb_ne = [];
ne_nb_ne = [];

d = dir('*-fracmb-pairs-resptype.mat');

for n = 1:length(d)

    filename = d(n).name;
    s = load(filename, 'rt');
    rt = s.rt;
    
    chan = [rt.chan];
    chan_unique = unique(chan);
    
    for i = 1:length(chan_unique)
        
        index = find( chan_unique(i) == chan );
        
        if length(index)>1
            
            cmb = nchoosek(index, 2); % determine all possible pairwise combinations
            
            [nr, nc] = size(cmb);
            
            for j = 1:nr
                
                index1 = cmb(j,1);
                index2 = cmb(j,2);
                
                position = [position; rt(index1).position];
                
                nb_nb_ne = [nb_nb_ne; rt(index1).nb_nb_ne_total rt(index2).nb_nb_ne_total];
                ne_nb_ne = [ne_nb_ne; rt(index1).ne_nb_ne_total rt(index2).ne_nb_ne_total];
                
            end % (for j)
            
        end % (if/else)
        
    end % (for i)

end % (for n)
return;



