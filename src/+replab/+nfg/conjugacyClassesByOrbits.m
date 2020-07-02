function classes = conjugacyClassesByOrbits(group)
% Generates the conjugacy classes of a group
%
% Note that it returns the results through the group nice monomorphism
%
% Based on the orbit computation algorithm from:
% Butler, Gregory. "Orbits and Schreier Vectors." FUNDAMENTAL ALGORITHMS FOR PERMUTATION GROUPS,
% by Gregory Butler, Springer-Verlag, 1991, pp. 57-59.
%
% Args:
%   group (`+replab.NiceFiniteGroup`): Nice finite group to compute the conjugacy classes of
%
% Returns:
%   classImages (cell(1,\*) of double(\*,\*) arrays): All conjugacy classes, as matrices
%                                                     with permutations in the class as rows

    I = replab.nfg.niceMonomorphismImages(group);
    ord = size(I, 1); % group order
    ds = size(I, 2); % domain size for images
    assert(ds < 65536, 'Domain size too big for naive enumeration');
    % we use a simple hash function that maps permutations to doubles
    % with a domain size < 2^16, that means that if h has values between -1023 and 1023,
    % the maximal value of the hash is 2^16*(2^16*2^10) = 2^36 which fits in a double
    h = randi([-1023 1023], ds, 1)-1;
    Ih = I * h;

    nG = group.nGenerators;
    gens = zeros(nG, ds);
    gensInv = zeros(nG, ds);
    for i = 1:nG
        gens(i,:) = group.niceMonomorphismImage(group.generator(i));
        gensInv(i,:) = group.niceMonomorphismImage(group.generatorInverse(i));
    end

    conjcl = zeros(1, ord);
    nConjCl = 0;
    for i = 1:ord
        if conjcl(i) == 0
            nConjCl = nConjCl + 1;
            conjcl(i) = nConjCl;
            toCheck = [i];
            while ~isempty(toCheck)
                g = I(toCheck(end), :);
                toCheck = toCheck(1:end-1);
                for j = 1:nG
                    gen = gens(j,:);
                    genInv = gensInv(j,:);
                    cj = genInv(g(gen)); % faster compose(genInv, g, gen)
                    f = find(Ih == cj*h);
                    if length(f) > 1
                        disp('bing');
                        % several rows have the same hash, so we look for an exact match
                        [~, loc] = ismember(cj, I(f,:), 'rows');
                        f = f(loc);
                    end
                    if conjcl(f) == 0
                        conjcl(f) = nConjCl;
                        toCheck = [toCheck f];
                    else
                        assert(conjcl(f) == nConjCl);
                    end
                end
            end
        end
    end
    classes = cell(1, nConjCl);
    for i = 1:nConjCl
        classes{i} = I(conjcl == i, :);
    end
end