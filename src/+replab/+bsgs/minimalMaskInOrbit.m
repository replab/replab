function mask = minimalMaskInOrbit(degree, generators, baseOrdering)
    if size(generators, 2) == 0
        mask = true(1, degree);
    else
        if nargin < 3
            baseOrdering = 1:degree;
        end
        mask = false(1, degree);
        done = false(1, degree);
        i = 1;
        while ~isempty(i)
            toCheck = i;
            minVal = baseOrdering(i);
            minInd = i;
            done(i) = true;
            while 1
                images = generators(toCheck, :);
                images = images(:);
                toCheck = images(~done(images));
                if isempty(toCheck)
                    break
                else
                    done(toCheck) = true;
                    [minVal1, minInd1] = min(baseOrdering(toCheck));
                    if minVal1 < minVal
                        minVal = minVal1;
                        minInd = toCheck(minInd1);
                    end
                end
            end
            mask(minInd) = true;
            i = find(~done, 1);
        end
        % simpler code: DEBUG
        orbits = replab.Partition.permutationsOrbits(generators').blocks;
        mask1 = false(1, degree);
        for i = 1:length(orbits)
            orbit = orbits{i};
            [~, ind] = min(baseOrdering(orbit));
            mask1(orbit(ind)) = true;
        end
        assert(isequal(mask, mask1));
    end
end
