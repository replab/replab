function [keyValuePairs exists oldValue] = keyValuePairsUpdate(keyValuePairs, key, newValue)
% Updates a sequence of key/value pairs
%
% Args:
%   keyValuePairs (cell(1,\*)): Key/value pairs (odd indices are charstring keys, even indices are values)
%   key (charstring): Key to add or update
%   newValue: Value to add or update
%
% Returns
% -------
%   keyValuePairs:
%     cell(1,\*): Updated list of key/value pairs
%   exists:
%     logical: Whether the key was already present
%   oldValue:
%     Previous value if already present
    [exists, ind] = ismember(key, keyValuePairs(1:2:end));
    if exists
        oldValue = keyValuePairs{ind*2};
        keyValuePairs{ind*2} = newValue;
    else
        keyValuePairs{1,end+1} = key;
        keyValuePairs{1,end+1} = newValue;
        oldValue = [];
    end
end
