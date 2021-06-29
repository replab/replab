function [q, u] = randomHouseholderMatrix(d, divisionRing)
% Returns a random Householder transformation over the real, complex, quaternion numbers
%
% Based on
% F. Mezzadri, "How to Generate Random Matrices from the Classical Compact Groups,"
% Notices of the AMS, vol. 54, no. 5, pp. 592–604, 2007.
%
% The Householder transformation is then given by matrix ``-conj(q)*(eye(d) - 2*u*u')``.
%
% Args:
%   d (integer): Dimension
%   divisionRing ('R', 'C', 'H'): Division ring
%
% Returns:
%   q: double or `+replab.H`
%     Phase factor
%   double(d,1) or `+replab.H`(d,1): Vector
    e1 = [1; zeros(d-1, 1)];
    switch divisionRing
      case 'R'
        v = randn(d, 1);
      case 'C'
        v = randn(d, 1) + 1i * randn(d, 1);
      case 'H'
        v = replab.H(randn(d, 1), randn(d, 1), randn(d, 1), randn(d, 1));
    end
    v = v / norm(v);
    x1 = abs(v(1));
    q = v(1)/x1; % or "sign(x1)"
    u = v + q * e1;
    u = u / norm(u);
end
