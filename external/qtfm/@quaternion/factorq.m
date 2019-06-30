function [f, g, theta] = factorq(q, mu, L)
% FACTORQ    Orthogonal factorization of a quaternion. The result is such
% that:
%
% f .* g = q, (factorization on the right), or
% g .* f = q, (factorization on the left),
%
% where in both cases g = exp(mu .* theta). Therefore g is a unit
% quaternion. The factor f will be orthogonal to the pure quaternion mu.

% Copyright (c) 2016 Stephen J. Sangwine and Todd A. Ell (This function was
% contributed to the QTFM toolbox based on mathematical results by Todd Ell
% 'Factoring Unit Quaternions', private communication, 11 September 2003.)
% See the file : Copyright.m for further details.

narginchk(3, 3), nargoutchk(0, 3)

if ~isa(q, 'quaternion') || ~isa(mu, 'quaternion')
    error('First two parameters must be quaternions')
end

if ~ispure(mu)
    error('Second parameter must be a pure quaternion.')
end

if ~ischar(L)
    error('Third parameter must be a character string');
end

mu = unit(mu); % Ensure mu has unit modulus.

if ispure(q)
    theta = pi ./ 2; % By definition if there is no scalar part.
else
    theta = atan2(scalar_product(mu, vector(q)), scalar(q));
end

g = exp(mu .* theta);

switch L
    case 'L'
        f = conj(g) .* q;
    case 'R'
        f = q .* conj(g);
    otherwise
        error('Third parameter must be ''L'' or ''R''')
end

end

% $Id: factorq.m 1020 2019-03-24 17:15:31Z sangwine $
