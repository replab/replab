function u = U(d)
% Returns the group of d x d unitary complex matrices
    u = replab.ClassicalCompactGroup(d, 'C', false);
end
