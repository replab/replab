function Q = QuaternionGroup
% Returns the quaternion group as a signed permutation groupq
%
% The quaternion group acts on the basis elements ``1``, ``j`, ``k``, ``l``.
%
% Returns:
%   `+replab.SignedPermutationGroup`: The quaternion group
    g1 = [-1 -2 -3 -4];
    gi = [2 -1 4 -3];
    gj = [3 -4 -1 2];
    Q = replab.SignedPermutationGroup.of(g1, gi, gj);
end
