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
    group = replab.CyclicGroup(n);
    classes = group.conjugacyClasses;
    % Check that conjugacy classes are sorted
    assert(isequal(cellfun(@(x) x.representative(1), classes), 1:n))
    chars = cell(n);
    for i = 1:n
        for j = 1:n
            chars{i, j} = sprintf('E(%d)^%d', n, mod((i-1) * (j-1), n));
        end
    end
    irrepExp = cellfun(@(x) sprintf('E(%d)^%d', n, x), num2cell(0:n-1), 'UniformOutput', false);
    ct = replab.CharacterTable.make(group, classes, [], chars, [], irrepExp);
end

