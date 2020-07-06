function rep = youngSeminormal(partition)
% Returns an irrep of S(n) according to the Young seminormal basis
%
% Follows the SAGE conventions
%
% http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/symmetric_group_representations.html
%
% Args:
%   partition (double): A Young diagram described by a row vector of row lengths
%
% Returns:
%   A `+replab.Rep` nonunitary irreducible representation
    n = sum(partition);
    Sn = replab.SymmetricGroup(n);
    [~, ~, NA, NB, ~, ~] = replab.sym.symIrrepImages(partition);
    d = size(NA, 1);
    if n == 1
        rep = Sn.repByImages('R', d, cell(1, 0), cell(1, 0));
        rep.isUnitary = true;
    elseif n == 2
        rep = Sn.repByImages('R', d, {NB}, {NB});
        rep.isUnitary = false;
    else
        NAinv = NA;
        for i = 1:n-2
            NAinv = NAinv * NA;
        end
        rep = Sn.repByImages('R', d, {NA NB}, {NAinv NB});
        rep.isUnitary = false;
    end
end
