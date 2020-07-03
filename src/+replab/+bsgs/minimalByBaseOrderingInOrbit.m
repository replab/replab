function mask = minimalByBaseOrderingInOrbit(degree, generators, baseOrdering)
    if size(generators, 2) == 0
        mask = true(1, degree);
    else
        orbits = replab.Partition.permutationsOrbits(generators').blocks;
        mask = false(1, degree);
        for i = 1:length(orbits)
            orbit = orbits{i};
            [~, ind] = min(baseOrdering(orbit));
            mask(orbit(ind)) = true;
        end
    end
end
