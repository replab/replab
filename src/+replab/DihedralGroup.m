function grp = DihedralGroup(n)
% Constructs the dihedral group of order ``2*n``
%
% This corresponds to the group of symmetries of the polygon with ``n`` vertices
%
% Args:
%   n (integer): Half the dihedral group order
%
% Returns:
%   `+replab.PermutationGroup`: The dihedral group permuting the vertices of the ``n``-gon
    if n <= 2
        grp = replab.SymmetricGroup(2);
    else
        g1 = fliplr(1:n);
        g2 = [2:n 1];
        grp = replab.PermutationGroup.of(g1, g2);
    end
end
