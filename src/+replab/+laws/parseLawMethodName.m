function [lawName sets] = parseLawMethodName(laws, methodName)
% Parses a law method name
%
% For the method name ``law_my_super_law_name_A1BC``, it will return
% ``lawName = 'my super law name'`` and ``domains = {laws.A1 laws.B laws.C}``.
%
% Args:
%   laws (`+replab.Laws`): Laws instance containing the method
%   methodName (charstring): Method name
%
% Returns
% -------
%   lawName:
%     charstring: The law name, camel case being transformed to space separated words
%   sets:
%     cell{1,:} of `+replab.Domain`: The domains corresponding to the type identifiers
    parts = strsplit(methodName, '_', 'CollapseDelimiters', 0);
    % we need at least the 'law' prefix, a 1-word law name, and a possibly empty type identifier string
    assert(length(parts) >= 3);
    assert(isequal(parts{1}, 'law'), 'Method name must start with ''law''');
    nameParts = parts(2:end-1);
    lawName = strjoin(nameParts, ' ');
    tis = parts{end}; % type identifier string
    check = regexp(tis, '^([A-Za-z][0-9]*)*$', 'once');
    types = regexp(tis, '([A-Za-z][0-9]*)', 'tokens');
    nTypes = length(types);
    sets = cell(1, nTypes);
    for i = 1:nTypes
        t = types{i};
        if isa(t, 'cell')
            t = t{1};
        end
        sets{i} = laws.(t);
    end
end
