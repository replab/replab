function res = isSimple(group)
% Computes whether the given group is simple
%
% Args:
%   `+replab.PermutationGroup`: Permutation group
%
% Returns:
%   logical: Whether the given group is simple
    if group.isTrivial
        res = false;
        return
    end
    C = group.conjugacyClasses;
    for i = 1:C.nClasses
        c = C.classes{i}.representative;
        if ~group.isIdentity(c)
            if group.normalClosure(group.subgroup({c})) ~= group
                res = false;
                return
            end
        end
    end
    res = true; % all conjugacy classes generate the full group
end