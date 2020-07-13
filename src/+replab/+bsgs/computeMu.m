function mu = computeMu(degree, l, base, g, resBasicOrbits, baseOrdering)
% Computes the ``mu`` bound in the backtracking tests
%
% See `.subgroupSearch`
    mu = degree + 2; % place holder for element < all others in base ordering
    for j = 1:l
        if any(base(l) == resBasicOrbits{j})
            candidate = g{j}(base(j));
            if baseOrdering(candidate) > mu
                mu = candidate;
            end
        end
    end
end
