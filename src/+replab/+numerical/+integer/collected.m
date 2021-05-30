function [elements, counts] = collected(list)
% For an unordered list of integers, return each integer with its number of occurences
%
% Args:
%   list (integer(1,\*)): List of integers
%
% Returns
% -------
%   elements: integer(1,\*)
%     Unique elements in list
%   counts: integer(1,\*)
%     Counts for each element
    list = list(:)';
    elements = unique(list);
    counts = sum(bsxfun(@eq, list', elements), 1);
end
