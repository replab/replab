function s = populateStruct(s, keyValuePairs)
% Populates a struct from key-value pairs
%
% Args:
%   s (struct): Struct with possible fields populated default values
%   keyValuePairs (cell(1,\*)): Additional arguments given as key-value pairs
%
% Returns:
%   struct: Struct with updated fields
    ind = 1;
    n = length(keyValuePairs);
    while ind <= n
        key = keyValuePairs{ind};
        assert(ischar(key), 'Named arguments must be key-value pairs, where the key in position %d is a charstring', ind);
        ind = ind + 1;
        assert(ind <= n, 'Named argument key %s must be followed by value');
        value = keyValuePairs{ind};
        ind = ind + 1;
        if ~isfield(s, key)
            error('No named argument %s exists', key);
        end
        s.(key) = value;
    end
end
