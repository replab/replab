function [s, rest] = populateStruct(s, keyValuePairs)
% Populates a struct from key-value pairs
%
% Raises:
%   An error if a key is not a charstring. An error if an unknown key is provided and less than two outputs
%   arguments are returned.
%
% Args:
%   s (struct): Struct with possible fields populated default values
%   keyValuePairs (cell(1,\*)): Additional arguments given as key-value pairs
%
% Returns
% -------
%   s:
%     struct: Struct with updated fields
%   rest:
%     cell(1,\*): Key-value pairs not part of the provided structure
    ind = 1;
    n = length(keyValuePairs);
    if nargout > 1
        rest = cell(1, 0);
    end
    while ind <= n
        key = keyValuePairs{ind};
        assert(ischar(key), 'Named arguments must be key-value pairs, where the key in position %d is a charstring', ind);
        ind = ind + 1;
        assert(ind <= n, 'Named argument key %s must be followed by value');
        value = keyValuePairs{ind};
        ind = ind + 1;
        if ~isfield(s, key)
            if nargout < 2
                error('No named argument %s exists', key);
            else
                rest{1, end+1} = key;
                rest{1, end+1} = value;
            end
        end
        s.(key) = value;
    end
end
