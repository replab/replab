function value = atlasEntries(newValue)
% Library of groups and character tables
%
% Args:
%   newValue (integer or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   cell(1,\*) of `+replab.AtlasEntry`: Stored atlas entries
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = cell(1, 0);
    end
    value = storedValue;
end
