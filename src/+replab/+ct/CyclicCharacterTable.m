function ct = CyclicCharacterTable(n)
% Generates the character table for the cyclic group Cn
%
% From Fässler, A., Stiefel, E., & Wong, B. D. (1992). Group theoretical methods and their applications.
% Boston: Birkhäuser, 22.
%
% Args:
%   n (integer): order of cyclic group
%
% Returns:
%   ct (`+replab.CharacterTable`)
    group = replab.PermutationGroup.cyclic(n);

    % Generate conjugacy class representatives in order gen^0, gen^1, gen^2, ...
    classreps = cell(1, n);
    gen = [2:n, 1];
    rep = 1:n;
    for i = 1:n
        classreps{i} = rep;
        rep = rep(gen);
    end
    classarray = cellfun(@(r) group.conjugacyClass(r), classreps, 'UniformOutput', false);
    classes = replab.ConjugacyClasses(group, classarray);

    % Generate irreps with images as increasing powers of E(n)
    w = replab.cyclotomic.E(n);
    irreps = arrayfun(@(x) group.repByImages('C', 1, 'images', {w^x}), 0:n-1, 'uniform', 0);

    % Generate characters with increasing powers of conjugacy classes and irreps
    chars = replab.cyclotomic.zeros(n, n);
    for i = 1:n
        for j = 1:n
            chars(i, j) = w^mod((i-1) * (j-1), n);
        end
    end

    ct = replab.CharacterTable(group, classes, chars, 'irreps', irreps);
end
