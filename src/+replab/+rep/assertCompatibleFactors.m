function b = assertCompatibleFactors(group, field, factors)
% Checks that the given representations are compatible in a direct sum or tensor product
%
% Args:
%   group (`+replab.CompactGroup`): Group on which the representations are defined
%   factors (cell(1,\*) of `+replab.Rep`): Representations of the given group
    assert(isa(group, 'replab.CompactGroup'));
    assert(ismember(field, {'R' 'C'}));
    assert((isempty(factors) || isrow(factors)) && iscell(factors));
    n = length(factors);
    for i = 1:n
        assert(isa(factors{i}, 'replab.Rep'));
        assert(group == factors{i}.group);
        assert(isequal(field, factors{i}.field));
    end
end
