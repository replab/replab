function rep = specht(partition)
    n = sum(partition);
    Sn = replab.Permutations(n);
    [SA, SB, ~, ~, ~, ~] = replab.nu.symIrrepImages(partition);
    d = size(SA, 1);
    if n == 2
        rep = replab.nu.RepByImages(Sn, 'R', d, SB, SB);
    else
        SAinv = SA;
        for i = 1:n-2
            SAinv = SAinv * SA;
        end
        rep = replab.nu.RepByImages(Sn, 'R', d, {SA SB}, {SAinv SB});
    end
end
