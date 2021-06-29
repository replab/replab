function grp = CyclicGroup(n)
% Describes cyclic permutations over n = "domainSize" elements

    warning('Deprecated. Use replab.PermutationGroup.cyclic(n) instead of replab.CyclicGroup(n)');
    grp = replab.PermutationGroup.cyclic(n);
end
