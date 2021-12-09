function b = isExact(mat)
% Returns whether the given matrix is exact: exact means either integer entries or cyclotomic
%
% Args:
%   mat (double(\*,\*) or `+replab.cyclotomic` (\*,\*)): Matrix
%
% Returns:
%   logical: Whether the matrix is exact
    b = isa(mat, 'replab.cyclotomic') || full(all(all(mat == round(mat))));
end
