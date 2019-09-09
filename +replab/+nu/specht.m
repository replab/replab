function rep = specht(partition)
% Returns an irrep of the symmetric group S(n) according to the Specht basis
%
% Follows the SAGE conventions
% 
% http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/symmetric_group_representations.html
%
% Args:
%   partition (double): A Young diagram described by a row vector of row lengths
%
% Returns:
%   A `replab.nu.Rep` nonunitary irreducible representation
    n = sum(partition);
    Sn = replab.Permutations(n);
    [SA, SB, ~, ~, ~, ~] = replab.nu.symIrrepImages(partition);
    d = size(SA, 1);
    if n == 1
        rep = replab.nu.RepByImages(Sn, 'R', d, {}, {});
    elseif n == 2
        rep = replab.nu.RepByImages(Sn, 'R', d, {SB}, {SB});
    else
        SAinv = SA;
        for i = 1:n-2
            SAinv = SAinv * SA;
        end
        rep = replab.nu.RepByImages(Sn, 'R', d, {SA SB}, {SAinv SB});
    end
end
