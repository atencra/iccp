function iccp_folder_plot_fra_sustained_resptype_raster

library('frabox');

[rtStr] = iccp_pairs_fra_resptype_files_to_struct;

ne = [rtStr.ne_nb_ne];
index = find(ne > 0.55);
neStr = rtStr(index);


for i = 1:length(neStr)

    fprintf('Processing %.0f of %.0f\n', i, length(neStr));

    file = neStr(i).file;
    indexNe = neStr(i).index;

    % Get FRA file and load the fra data structure
    index = findstr(file, '-resptype.mat');
    fraFile = sprintf('%s.mat', file(1:index-1));
    load(fraFile);

    % Plot the raster for the FRA
    plot_freq_resp_area_raster(fra(indexNe));

end % (for i)






function [rtStr] = iccp_pairs_fra_resptype_files_to_struct

position = [];
nb_nb_ne = [];
ne_nb_ne = [];

d = dir('*-fracmb-pairs-resptype.mat');

rtStr = [];

for n = 1:length(d)

    filename = d(n).name;
    s = load(filename, 'rt');
    rt = s.rt;

    for j = 1:length(rt)

        rtStrTemp.file = filename;
        rtStrTemp.index = j;
        rtStrTemp.nb_nb_ne = rt(j).nb_nb_ne_total;
        rtStrTemp.ne_nb_ne = rt(j).ne_nb_ne_total;

        rtStr = [rtStr rtStrTemp];
    end % (for j)

end % (for n)

return;


