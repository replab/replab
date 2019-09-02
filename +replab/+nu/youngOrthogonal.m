function rep = youngOrthogonal(partition)
    n = sum(partition);
    Sn = replab.Permutations(n);
    [~, ~, ~, ~, OA, OB] = replab.nu.symIrrepImages(partition);
    d = size(OA, 1);
    if n == 2
        rep = Sn.rep('R', d, OB);
    else
        rep = Sn.rep('R', d, {OA OB});
    end
end
