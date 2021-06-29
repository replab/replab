function o = SO(d)
% Returns the group of d x d orthonormal matrices with determinant one
    o = replab.ClassicalCompactGroup(d, 'R', true);
end
