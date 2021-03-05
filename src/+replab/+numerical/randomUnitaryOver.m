function [U, detU] = randomUnitaryOver(d, divisionRing)
% Samples an orthogonal/unitary/compact symplectic according to the Haar measure
%
% Based on
% F. Mezzadri, "How to Generate Random Matrices from the Classical Compact Groups,"
% Notices of the AMS, vol. 54, no. 5, pp. 592–604, 2007.
%
% See also `<http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf>_`
%
% An alternative algorithm would be given by:
% `<http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf>_`
%
% Args:
%   d (integer): Dimension
%   divisionRing ('R', 'C', 'H'): Division ring for the matrix coefficients
%
% Returns
% -------
%   U: double(d,d) or `+replab.H`(d,d)
%     Random matrix
%   detU: double
%     Matrix determinant (or ``[]`` if divisionRing = 'H')
    [q, u] = replab.numerical.randomHouseholderMatrix(1, divisionRing);
    % The first step is simplified
    % -conj(q)*(eye(d) - 2*u*u') = -conj(q)*(1 - 2) = conj(q)
    U = conj(q);
    if divisionRing == 'H'
        detU = [];
    else
        detU = conj(q);
    end
    % Now we iterate over the cosets of subgroups
    for i = 2:d
        [q, u] = replab.numerical.randomHouseholderMatrix(i, divisionRing);
        % Below, commented, the original computation, which we changed for speed reasons
        % U = [1, zeros(1, i-1); zeros(i-1, 1), U];
        % U = -conj(q)*(U - 2*u*(u'*U));
        % or U = -conj(q)*(eye(d) - 2*u*u')*U
        uU = u(2:end)'*U;
        if divisionRing ~= 'H'
            % we keep track and update the determinant, as the reflection part (1 - 2*u*u') has det -1
            detU = detU * conj(q)^i * (-1)^(i-1);
        end
        U = -conj(q)*[1-2*u(1)*conj(u(1)), -2*u(1)*uU; -2*u(2:end)*conj(u(1)), U - 2*u(2:end)*uU];
    end
    % One step of Newton iteration to force unitary (helps!)
    N = U'*U;
    P = U*N/2;
    U = 2*U + P*N - 3*P;
end
