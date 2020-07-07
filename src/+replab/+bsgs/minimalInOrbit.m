function minInd = minimalInOrbit(degree, generators, b, baseOrdering)
    if size(generators, 2) == 0
        minInd = b;
    else
        if nargin < 4
            baseOrdering = 1:degree;
        end
        % simpler code:
        % orbits = replab.Partition.permutationsOrbits(generators');
        % m = min(orbits.blocks{orbits.blockIndex(b)})
        done = false(1, degree);
        toCheck = b;
        minVal = baseOrdering(b);
        minInd = b;
        done(b) = true;
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
                    minInd = toCheck(minInd1);
                end
            end
        end
    end
end
