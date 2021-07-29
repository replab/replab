function grp = D(n)
% Constructs the dihedral group of order ``2*n``
%
% This corresponds to the group of symmetries of the polygon with ``n`` vertices
%
% Args:
%   n (integer): Half the dihedral group order
%
% Returns:
%   `+replab.PermutationGroup`: The dihedral group permuting the vertices of the ``n``-gon
    grp = replab.PermutationGroup.dihedral(n);
end
