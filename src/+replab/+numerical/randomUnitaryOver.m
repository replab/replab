function U = randomUnitaryOver(d, divisionRing)
% Samples an orthogonal/unitary/compact symplectic according to the Haar measure
%
% Based on
% F. Mezzadri, "How to Generate Random Matrices from the Classical Compact Groups,"
% Notices of the AMS, vol. 54, no. 5, pp. 592–604, 2007.
%
% An alternative algorithm would be given by:
% `<http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf>_`
%
% Args:
%   d (integer): Dimension
%   divisionRing ('R', 'C', 'H'): Division ring for the matrix coefficients
%
% Returns:
%   double(d,d) or `+replab.H`(d,d): Random matrix
    [q, u] = replab.numerical.randomHouseholderMatrix(1, divisionRing);
    U = conj(q); % -conj(q)*(eye(d) - 2*u*u') = -conj(q)*(1 - 2) = conj(q)
    for i = 2:d
        [q, u] = replab.numerical.randomHouseholderMatrix(i, divisionRing);
        %U = [1, zeros(1, i-1); zeros(i-1, 1), U];
        %U = -conj(q)*(U - 2*u*(u'*U));
        %or U = -conj(q)*(eye(d) - 2*u*u')*U
        uU = u(2:end)'*U;
        U = -conj(q)*[1-2*u(1)*conj(u(1)), -2*u(1)*uU; -2*u(2:end)*conj(u(1)), U - 2*u(2:end)*uU];
    end
    % One step of Newton iteration to force unitary (helps!)
    N = U'*U;
    P = U*N/2;
    U = 2*U + P*N - 3*P;
end
