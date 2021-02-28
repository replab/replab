function u = SU(d)
% Returns the group of d x d unitary complex matrices with determinant one
    u = replab.ClassicalCompactGroup(d, 'C', true);
end
