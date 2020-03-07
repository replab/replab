function c = fusionOptional(a, b)
% Merges optional values, checking for consistency
%
% - If both ``a`` and ``b`` are empty, returns ``[]``.
% - If only one of them is empty, return the other.
% - If both ``a`` and ``b`` are non-empty and equal, return that value.
% - If both ``a`` and ``b`` are non-empty and differ, raises an error.
%
% Raises:
%   An error if ``a`` and ``b`` are non-empty and differ.
%
% Args:
%   a: First value
%   b: Second value
%
% Returns:
%   A non-empty value if ``a`` and ``b`` are consistent, otherwise ``[]``
    if isempty(a)
        if isempty(b)
            c = [];
        else
            c = b;
        end
    else
        if isempty(b)
            c = a;
        else
            assert(isequal(a, b));
            c = a;
        end
    end
end
