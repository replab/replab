function minInd = minimalInOrbit(degree, generators, b, ordering)
% Returns the minimal point in the orbit of a given point according to a custom ordering
%
% The comparison function ``x <bo y`` is given by ``x <bo y = ordering(x) < ordering(y)``.
%
% Args:
%   generators (integer(\*,\*)): Group generators given as matrix columns
%   b (integer): Point to compute the orbit of
%   ordering (integer(1,\*), optional): Ordering, default value ``[1 2 ... n]``
%
% Returns:
%   integer: Minimal point in the orbit
    if size(generators, 2) == 0
        minInd = b;
    else
        if nargin < 4
            ordering = 1:degree;
        end
        % simpler code:
        % orbits = replab.Partition.permutationsOrbits(generators');
        % m = min(orbits.blocks{orbits.blockIndex(b)})
        done = false(1, degree);
        toCheck = b;
        minVal = ordering(b);
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
                [minVal1, minInd1] = min(ordering(toCheck));
                if minVal1 < minVal
                    minInd = toCheck(minInd1);
                end
            end
        end
    end
end
