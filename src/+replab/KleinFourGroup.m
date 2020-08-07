function grp = KleinFourGroup()
% Constructs the Klein Four-Group
%
% This corresponds to the symmetry group of a non-square rectangle
%
% Returns:
%   `+replab.KleinFourGroup`: The Klein Four-Group
    g1 = [2,1,4,3];
    g2 = [3,4,1,2];
    grp = replab.PermutationGroup.of(g1, g2);
end
    