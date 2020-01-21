function rep = youngOrthogonal(partition)
% Returns an irrep of S(n) according to the Young orthogonal basis
%
% Follows the SAGE conventions
%
% http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/symmetric_group_representations.html
%
% Args:
%   partition (double): A Young diagram described by a row vector of row lengths
%
% Returns:
%   A `+replab.Rep` orthonormal irreducible representation
    n = sum(partition);
    Sn = replab.Permutations(n);
    [~, ~, ~, ~, OA, OB] = replab.sym.symIrrepImages(partition);
    d = size(OA, 1);
    if n == 1
        rep = Sn.repByImages('R', d, true, {}, {});
    elseif n == 2
        rep = Sn.repByImages('R', d, true, {OB});
    else
        rep = Sn.repByImages('R', d, true, {OA OB});
    end
end
