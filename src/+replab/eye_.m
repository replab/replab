function I = eye_(d)
% Returns the identity matrix of a specified dimension 
%
% Uses `replab.Settings.sparseCutoff` to decide whether the
% returned matrix is sparse.
%
% Args:
%   d (integer): Dimension of the matrix
%
% Returns:
%   double matrix: A `d` by `d` identity matrix
    if d >= replab.Settings.sparseCutoff
        I = speye(d);
    else
        I = eye(d);
    end
end
