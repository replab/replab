function classes = conjugacyClassesByOrbits(group)
% Generates the conjugacy classes of a group
%
% Based on the orbit computation algorithm from:
% Butler, Gregory. "Orbits and Schreier Vectors." FUNDAMENTAL ALGORITHMS FOR PERMUTATION GROUPS,
% by Gregory Butler, Springer-Verlag, 1991, pp. 57-59.
%
% Args:
%   group (`+replab.PermutationGroup`): Permutation group to compute the conjugacy classes of
%
% Returns:
%   classImages (cell(1,\*) of double(\*,\*) arrays): All conjugacy classes, as matrices
%                                                     with permutations in the class as columns
    I = group.chain.allElements;
    ds = size(I, 1); % domain size for images
    ord = size(I, 2); % group order
    assert(ds < 65536, 'Domain size too big for naive enumeration');
    % we use a simple hash function that maps permutations to doubles
    % with a domain size < 2^16, that means that if h has values between -2^20+1 and 2^20-1,
    % the maximal value of the hash is 2^16*(2^16*2^20) = 2^52 which fits in a double
    h = randi([-2^20+1 2^20-1], 1, ds);
    Ih = h * I;

    % sorts the hash values and the matrix of group elements
    % so that we can perform a binary search later
    [~, ind] = sort(Ih);
    I = I(:, ind);
    Ih = Ih(:, ind);

    % tests the existence of the binary search built-in
    has_ismembc2 = exist('ismembc2') > 0;
    nG = group.nGenerators;
    gens = zeros(ds, nG);
    gensInv = zeros(ds, nG);
    for i = 1:nG
        gens(:,i) = group.generator(i);
        gensInv(:,i) = group.generatorInverse(i);
    end

    conjcl = zeros(1, ord);
    nConjCl = 0;
    for i = 1:ord
        if conjcl(i) == 0
            nConjCl = nConjCl + 1;
            conjcl(i) = nConjCl;
            toCheck = [i];
            while ~isempty(toCheck)
                g = I(:, toCheck(end));
                toCheck = toCheck(1:end-1);
                for j = 1:nG
                    gen = gens(:,j);
                    genInv = gensInv(:,j);
                    cj = genInv(g(gen)); % faster compose(genInv, g, gen)
                                         % f = ismembc2(h*cj, Ih);
                    hval = h*cj;
                    if has_ismembc2
                        % perform binary search; this function returns the last element
                        % which matches
                        last = ismembc2(hval, Ih);
                        % then we need to find if other elements before match as well
                        before = last - 1;
                        while before > 0 && Ih(before) == hval
                            before = before - 1;
                        end
                        f = before+1:last;
                    else
                        % fallback on the default find
                        f = find(Ih == hval);
                    end

                    if length(f) > 1
                        % several rows have the same hash, so we look for an exact match
                        loc = replab.util.findRowInMatrix(cj', I(:,f)');
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
        classes{i} = I(:, conjcl == i);
    end
end
