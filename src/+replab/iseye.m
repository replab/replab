function b = iseye(M)
% Tests whether a given matrix is exactly the identity
%
% Args:
%   M (double matrix, can be sparse): Arbitrary matrix
%
% Returns:
%   logical: True if the matrix is exactly the identity matrix, false otherwise
    if size(M, 1) == size(M, 2) && isdiag(M)
        b = all(diag(M) == 1);
    else
        b = false;
    end
end
