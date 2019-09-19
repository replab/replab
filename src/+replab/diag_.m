function D = diag_(vec)
% Returns the diagonal matrix with given elements on the diagonal
%
% Implements a subset of the functionality of the `diag` Octave/Matlab function.
%
% Makes sure to use sparse matrices starting from `replab.Settings.sparseCutoff`.
% Args:
%   vec (row vector of double): Diagonal elements
%
% Returns:
%   double matrix: Diagonal matrix
    d = length(vec);
    if d >= replab.Settings.sparseCutoff
        D = sparse(1:d, 1:d, vec);
    else
        D = diag(vec);
    end
end
