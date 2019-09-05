function [s1 s2] = size(self, dim)
    % returns the matrix size
    if (nargin >= 2) && (dim ~= 1) && (dim ~= 2)
        error('Wrong dimension in gem::size');
    end

    s1 = self.dim*[1 1];

    if nargin >= 2
        s1 = s1(dim);
    elseif nargout == 2
        s2 = s1(2);
        s1 = s1(1);
    end
end
