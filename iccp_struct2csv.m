function iccp_struct2csv(data,csvfile)
% iccp_struct2csv Write struct data to csv file 
% 
%
%    data : struct holding paired data. Each field of data holds an
%    Nx1 vector or an Nx2 array.
%
%    data is obtained from functions like:
%
%       data = iccp_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
%       data = iccp_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs)
%       data = iccp_plot_strf_bmf_crosscorr_pos_neg(ccpairs)
%       data = iccp_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs)
%
%
%    Data is written to a CSV file, called data.csv, with one row per pair 
%    of neurons.
%
%    Headers are supplied on the first line of the csv file. Headers are 
%    the fields of the struct. For the Nx2 arrays, each column is written 
%    to a column of the csv, and the headers have a '1' or a '2' appended 
%    to the name.
%
%    csvfile = iccp_strf_si_ccc_to_csv(data,csvfile) saves the data in the
%    file csvfile.csv.


narginchk(1,2);

if ( nargin == 1 )
    csvfile = 'data.csv';
end

if ( nargin == 2 )
    if ( isempty(findstr(csvfile,'.csv')) )
        csvfile = sprintf('%s.csv',csvfile);
    end
end

fields = fieldnames(data);


fid = fopen(csvfile,'w');

for i = 1:length(fields)

    d = eval(sprintf('data.%s',fields{i}));
    [~,nc] = size(d);

    if ( nc == 1 )
        if ( i < length(fields) )
            fprintf(fid, '%s,',fields{i});
        else
            fprintf(fid, '%s\n',fields{i});
        end

    else
        if ( i < length(fields) )
            for j = 1:nc
                fprintf(fid, '%s%.0f,',fields{i},j);
            end
        else
            for j = 1:nc
                if ( j < nc )
                    fprintf(fid, '%s%.0f,',fields{i},j);
                else
                    fprintf(fid, '%s%.0f\n',fields{i},j);
                end
            end % (for j)
        end

    end

end % (for i)



% Now write out each row to the csv file
dataMat = [];
for i = 1:length(fields)
    d = eval(sprintf('data.%s',fields{i}));
    dataMat = [dataMat d];
end % (for i)



[nr,nc] = size(dataMat);

for i = 1:nr

    row = dataMat(i,:);

    for j = 1:nc
        if ( j < nc )
            fprintf(fid, '%.5f,', row(j));
        else
            fprintf(fid, '%.5f\n', row(j));
        end
    end

end % (for i)


fclose(fid);

   
return;




