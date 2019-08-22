function I = decomposeSymmetric(rep)
% Decomposes a representation of the symmetric group
%
% Args:
%   rep: Representation of the symmetric group S(n).
%        rep.group must of type 'replab.PermutationGroup', with ``rep.group.domainSize == n``
%        and ``rep.group.order == factorial(vpi(n))``
%
% Returns:
%   Either an ``replab.Irreducible`` instance, or ``[]`` if the representation
%   could not be decomposed.
    replab.Dispatch.assert(isa(rep.group, 'replab.PermutationGroup'), 'Not a permutation group');
    replab.Dispatch.assert(rep.overR, 'Not a real representation');
    n = rep.group.domainSize;
    d = rep.dimension;
    assert(rep.group.order == factorial(vpi(n)));
    % the two generators of the symmetric group we use in our canonical presentation
    a = [2:n 1];
    b = [2 1 3:n];
    try
        % find permutation representation images, throw exception if the rep is not a permutation rep
        a1 = replab.Permutations.fromMatrix(rep.image(a));
        b1 = replab.Permutations.fromMatrix(rep.image(b));
    catch
        replab.Dispatch.error('Not a permutation representation');
    end
    % find the transitive orbits
    P = replab.Partition.permutationsOrbits([a1; b1]);
    blockSizes = cellfun(@length, P.blocks);
    orderedToUnordered = [P.blocks{:}];
    unorderedToOrdered = replab.Permutations(length(orderedToUnordered)).inverse(orderedToUnordered);
    shift = 0;
    S = struct;
    Y = struct;
    for i = 1:P.nBlocks
        ordered = shift+(1:blockSizes(i));
        ai = unorderedToOrdered(a1(orderedToUnordered(ordered)))-shift;
        bi = unorderedToOrdered(b1(orderedToUnordered(ordered)))-shift;
        [Ui Yi] = replab.rep.decomposeSymmetricTransitive(n, ai, bi);
        for j = 1:length(Ui)
            Uij = Ui{j};
            Yij = Yi{j};
            dij = size(Uij, 1);
            Uij1 = zeros(dij, d);
            Uij1(:, orderedToUnordered(ordered)) = Uij;
            sub = replab.Irrep(rep, Uij1, replab.DivisionAlgebra.real);
            key = Yij.asIdentifier;
            if isfield(S, key)
                copies = S.(key);
            else
                copies = {};
                Y.(key) = Yij;
            end
            copies{1, end+1} = sub;
            S.(key) = copies;
        end
        shift = shift + blockSizes(i);
    end
    keys = fieldnames(Y);
    youngs = struct2cell(Y);
    isotypics = struct2cell(S);
    % reverse lexicographic sort
    p = replab.Permutations.sorting(youngs, @(x, y) replab.rep.lexicographicCompare(x.partition, y.partition) == -1);
    isotypics = cellfun(@(subs) replab.Isotypic(rep, subs), isotypics(p), 'uniform', 0);
    I = replab.Irreducible(rep, isotypics);
end
