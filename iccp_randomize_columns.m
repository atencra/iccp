function [absc, ord] = iccp_randomize_columns(a,b)
% iccp_randomize_columns Randomize pairwise data for plots
% 
%     [absc, ord] = iccp_randomize_columns(a,b)
% 
%     For data vectors a and b, where a and b have the same number of elements,
% 
%     [absc, ord] = iccp_randomize_columns(a,b) randomly assigns the value
%     in a(i,1) to either absc or ord. It then assigns b(i,2) to ord or abs.
% 
%     absc, ord: same length as a,b, but with values randomly chosen from each
%     row of a and b.



if nargin == 1
    [nr, nc] = size(a);
    if nc == 2
        b = a(:,2);
        a = a(:,1);
    else
        error('Wrong input size.');
    end
end
a = a(:);
b = b(:);


absc = zeros(size(a));
ord = zeros(size(b));

for i = 1:length(absc)
    r = rand(1) - 0.5;
    if r < 0
        absc(i) = a(i);
        ord(i) = b(i);
    else
        absc(i) = b(i);
        ord(i) = a(i);
    end
end



return
