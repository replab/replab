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
%   A `replab.Rep` orthonormal irreducible representation
    n = sum(partition);
    Sn = replab.Permutations(n);
    [~, ~, ~, ~, OA, OB] = replab.nu.symIrrepImages(partition);
    d = size(OA, 1);
    if n == 1
        rep = replab.nu.RepByImages(Sn, 'R', d, {}, {});
    elseif n == 2
        rep = Sn.rep('R', d, {OB});
    else
        rep = Sn.rep('R', d, {OA OB});
    end
end
