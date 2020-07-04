function nu = computeNu(degree, l, sortedOrbits, orbits, resBasicOrbits)
    idx = length(orbits{l}) + 2 - length(resBasicOrbits{l});
    if idx > length(sortedOrbits{l})
        nu = degree + 1; % place holder for element > all others in base ordering
    else
        nu = sortedOrbits{l}(idx);
    end
end
