function ct = CyclicCharacterTable(n)
% Generates the character table for the cyclic group Cn
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
    base = ['E(', num2str(n), ')'];
    chars = cell(n);
    for i = 1:n
        for j = 1:n
            chars{i, j} = [base, '^', num2str(mod((i-1) * (j-1), n))];
        end
    end
    irrepExp = cellfun(@(x) [base, '^', num2str(x)], num2cell(0:n-1), 'UniformOutput', false);
    ct = replab.CharacterTable.make(group, classes, [], chars, [], irrepExp);
end

