function s = S(n)
% Returns the symmetric group acting on n elements
%
% Alias for replab.Permutations(n)
%
% >> replab.S(3).order
% ans =
%     6
    s = replab.Permutations(n);
end
