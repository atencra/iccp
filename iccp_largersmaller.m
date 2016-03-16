function [larger, smaller] = iccp_largersmaller(a,b)

      larger = zeros( size(a) );
      smaller = zeros(size(b) );

      for ii = 1:length(a)
         if ( a(ii) > b(ii) )
            larger(ii) = a(ii);
            smaller(ii) = b(ii);
         else
            larger(ii) = b(ii);
            smaller(ii) = a(ii);
         end
      end % (for i)
return
