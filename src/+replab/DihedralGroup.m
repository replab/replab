function grp = DihedralGroup(n)
    if n <= 2
        grp = replab.SymmetricGroup(2);
    else
        g1 = [2:n 1];
        g2 = fliplr(1:n);
        grp = replab.PermutationGroup.of(g1, g2);
    end
end
