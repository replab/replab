function [A, B] = pairOfIntegerInverses(d, ct)
% Returns a pair of integer matrices that are inverse of each other
%
% Args:
%   d (integer): Dimension of the matrix
%   ct (integer): Constant integer factor parameterizing invertible matrices
%
% Returns
% -------
%   A: integer(d,d)
%     First matrix
%   B: integer(d,d)
%     Second matrix such that ``A * B == eye(d)``
    A = zeros(d,d);
    B = zeros(d,d);
    for i = 1:d
        for j = 1:d
            A(i,j) = mynchoosek(ct+j-1, i-1);
            for l = 0:d-j
                B(i,j) = B(i,j)+(-1)^(i+j)*mynchoosek(ct+l-1,l)*mynchoosek(l+j-1,i-1);
            end
        end
    end
end

function c = mynchoosek(p, q)
    if p == -1 && q == 0
        c = 1;
    elseif p < q || q < 0
        c = 0;
    else
        c = nchoosek(p, q);
    end
end
