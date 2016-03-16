function [rpop_med, rpop_ci, pval] = iccp_pairwise_corr_rand_test(data1, data2, logtransform)

if ( nargin == 2 )
    logtransform = 0;
end

rpop = zeros(1,1000);
for i = 1:1000
    [absc, ord] = iccp_randomize_columns(data1, data2);

    if ( logtransform )
        x = log10(absc);
        y = log10(ord);
        index = ~isnan(x) & ~isnan(y);
        x = x(index);
        y = y(index);
    else
        index = ~isnan(absc) & ~isnan(ord);
        x = absc(index);
        y = ord(index);
    end

    [r,p] = corrcoef(x, y);
    rpop(i) = r(2);
end % (for i)


n = length(data1);
data = [data1(:); data2(:)];
index_total = 1:length(data);

rpop_rand = zeros(1,1000);
for i = 1:1000
    index_absc = randperm(length(data), n);
    index_ord = setdiff(index_total, index_absc);
    absc = data(index_absc);
    ord = data(index_ord);

    if ( logtransform )
        x = log10(absc);
        y = log10(ord);
        index = ~isnan(x) & ~isnan(y);
        x = x(index);
        y = y(index);
    else
        index = ~isnan(absc) & ~isnan(ord);
        x = absc(index);
        y = ord(index);
    end

    [r,p] = corrcoef(x, y);
    rpop_rand(i) = r(2);
end % (for i)


rpop_sort = sort(rpop);
rpop_ci(1) = rpop_sort(25);
rpop_ci(2) = rpop_sort(975);
rpop_med = median(rpop);


index = find(rpop_rand > rpop_med);
if isempty(index)
    pval = 0;
else
    pval = length(index) / length(rpop_rand);
end


% figure;
% subplot(2,1,1);
% hist(rpop);
% subplot(2,1,2);
% hist(rpop_rand)

return;



