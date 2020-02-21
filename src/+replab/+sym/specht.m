function rep = specht(partition)
% Returns an irrep of the symmetric group S(n) according to the Specht basis
%
% Follows the SAGE conventions
%
% `<http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/symmetric_group_representations.html>`_
%
% Args:
%   partition (double): A Young diagram described by a row vector of row lengths
%
% Returns:
%   A `+replab.Rep` nonunitary irreducible representation
    n = sum(partition);
    Sn = replab.Permutations(n);
    [SA, SB, ~, ~, ~, ~] = replab.sym.symIrrepImages(partition);
    d = size(SA, 1);
    if n == 1
        rep = Sn.repByImages('R', d, {}, {});
    elseif n == 2
        rep = Sn.repByImages('R', d, {SB}, {SB});
    else
        SAinv = SA;
        for i = 1:n-2
            SAinv = SAinv * SA;
        end
        rep = Sn.repByImages('R', d, {SA SB}, {SAinv SB});
    end
    rep.isUnitary = false;
    rep.isIrreducible = true;
end
