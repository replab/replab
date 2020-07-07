function grp = DihedralGroup(n)
    if n <= 2
        grp = replab.SymmetricGroup(2);
    else
        g1 = fliplr(1:n);
        g2 = [2:n 1];
        grp = replab.PermutationGroup.of(g1, g2);
    end
end
