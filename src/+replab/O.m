function o = O(d)
% Returns the group of d x d orthonormal matrices
    o = replab.ClassicalCompactGroup(d, 'R', false);
end
