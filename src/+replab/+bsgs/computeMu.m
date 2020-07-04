function mu = computeMu(degree, l, base, g, resBasicOrbits, baseOrdering)
    mu = degree + 2; % place holder for element < all others in base ordering
    for j = 1:l
        if ismember(base(l), resBasicOrbits{j})
            candidate = g{j}(base(j));
            if baseOrdering(candidate) > mu
                mu = candidate;
            end
        end
    end
end
