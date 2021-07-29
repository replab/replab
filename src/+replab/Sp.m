function G = Sp(d)
% Returns the compact symplectic group whose elements are unitary quaternion matrices
    G = replab.ClassicalCompactGroup(d, 'H', false);
end
