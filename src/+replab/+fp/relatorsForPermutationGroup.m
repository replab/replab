function relators = relatorsForPermutationGroup(group)
% Computes relators for a finite group
%
% Note: calls the GAP system internally
%
% Args:
%   group (`+replab.PermutationGroup`): Permutation group to find the relators of
%
% Returns:
%   cell(1,\*) of integer(1,\*): Relators given as explicit words
    relators = cell(1, 0);
    if group.isTrivial
        return
    end
    nG = group.nGenerators;
    grp_s = group.subgroup({group.generator(1)});
    sub_s = group.trivialSubgroup;
    s = 1;
    while true
        relators = replab.fp.CosetTable.presentation(grp_s, sub_s, relators);
        s = s + 1;
        if s > nG
            break
        end
        sub_s = grp_s;
        grp_s = group.subgroup(group.generators(1:s));
    end
end
