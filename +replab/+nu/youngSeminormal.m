function rep = youngSeminormal(partition)
    n = sum(partition);
    Sn = replab.Permutations(n);
    [~, ~, NA, NB, ~, ~] = replab.nu.symIrrepImages(partition);
    d = size(NA, 1);
    if n == 2
        rep = replab.nu.RepByImages(Sn, 'R', d, NB, NB);
    else
        NAinv = NA;
        for i = 1:n-2
            NAinv = NAinv * NA;
        end
        rep = replab.nu.RepByImages(Sn, 'R', d, {NA NB}, {NAinv NB});
    end
end
