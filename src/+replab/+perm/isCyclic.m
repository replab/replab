function res = isCyclic(group)
% Computes whether the given group is cyclic
%
% Args:
%   `+replab.PermutationGroup`: Permutation group
%
% Returns:
%   logical: Whether the given group is cyclic
    if group.nGenerators <= 1
        res = true;
    elseif ~group.isCommutative
        res = false;
    else
        pds = factor(group.order);
        if length(pds) == 1 % group of prime order is cyclic
            res = true;
            return
        end
        assert(all(pds <= 2^53-1)); % to be sure (otherwise can a BSGS have been computed?)
        pds = unique(double(pds));
        pds = double(pds);
        for p = pds
            newGens = cellfun(@(g) group.composeN(g, p), group.generators, 'uniform', 0);
            newGens = newGens(~cellfun(@(g) group.isIdentity(g), newGens));
            if group.subgroup(newGens).order*p ~= group.order
                res = false;
                return
            end
        end
        res = true;
    end
end
