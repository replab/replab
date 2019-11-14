classdef FullSymmetryShape < replab.tensor.Shape
    
   methods
       
       function orbit = subOrbit(self, sub)
           [u, ~, J] = unique(sub);
           % Code adapted from https://www.mathworks.com/matlabcentral/newsreader/view_thread/164470
           p = u(up(J, length(sub)));
           function p = up(J, n)
               ktab = histc(J,1:max(J));
               l = n;
               p = zeros(1, n);
               s = 1;
               for i = 1:length(ktab)
                   k = ktab(i);
                   c = nchoosek(1:l, k);
                   m = size(c,1);
                   [t, ~] = find(~p.');
                   t = reshape(t, [], s);
                   c = t(c,:)';
                   s = s*m;
                   r = repmat((1:s)',[1 k]);
                   q = accumarray([r(:) c(:)], i, [s n]);
                   p = repmat(p, [m 1]) + q;
                   l = l - k;
               end
           end
       end
       
   end
   
   methods (Static)
      
       function dim = computeDimension(d, n)
            dim = nchoosek(n + d - 1, n);
       end
           
   end
   
end
