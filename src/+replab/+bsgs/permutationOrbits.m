function [orbits orbitSizes orbitIndex] = permutationOrbits(g)
    n = length(g);
    done = false(1, n);
    orbits = cell(1, 0);
    orbitSizes = zeros(1, n);
    orbitIndex = zeros(1, n);
    ind = 1;
    while ~all(done)
        b0 = find(~done, 1);
        orbit = b0;
        b = g(b0);
        while b ~= b0
            orbit = [orbit b];
            b = g(b);
        end
        orbitIndex(orbit) = ind;
        ind = ind + 1;
        orbitSizes(orbit) = length(orbit);
        done(orbit) = true;
        orbits{1,end+1} = orbit;
    end
end
