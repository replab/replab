function [U Y] = decomposeSymmetricTransitive(n, a, b)
% Decomposes a transitive permutation representation of the symmetric group
%
% Args:
%   n: Domain size of the symmetric group S(n), so that the order is n!
%   a: Permutation image of the cycle generator [2:n 1]
%   b: Permutation image of the transposition generator [2 1 3:n]
%
% Returns
% -------
%   U:
%     cell array of basis matrices, on matrix per irrep, with row basis vectors
%   Y:
%     row cell array of Young diagrams corresponding to the identified irreps
    if length(a) == 1 && length(b) == 1
        % trivial representation
        U = {1};
        Y = {replab.young.Diagram(n)};
        return
    end
    switch n
      case 2
        assert(length(a) == 2 && length(b) == 2);
        U = {[1  1] [1 -1]};
        Y = {replab.young.Diagram(2) replab.young.Diagram([1 1])};
      otherwise
        error('replab:Dispatch:tryNext');
    end
end
