function grp = CyclicGroup(n)
    warning('Deprecated. Use replab.PermutationGroup.cyclic(n) instead of replab.CyclicGroup(n)');
    grp = replab.PermutationGroup.cyclic(n);
end
