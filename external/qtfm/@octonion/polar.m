function [t, h, n] = polar(q, L)
% POLAR decomposition of a bioctonion. Computes a factorisation such that
% t .* h .* n = q (right polar decomposition, L == 'R' or omitted), or
% h .* t .* n = q (left polar decomposition, L == 'L').
% n is the semi-norm of q, and if omitted will trigger a warning if q is
% not of unit modulus. n is complex, in general.

% t is the trigonometric factor (a real octonion), invariant to the
% ordering. It will have a real axis and angle.
% h is the hyperbolic factor (a bioctonion), dependent on the ordering.
% It will have an imaginary axis and real angle. NB The axis and angle
% functions will return a real axis and imaginary angle (since the complex
% root of -1 commutes with the axis, it can be attached to either without
% invalidating the result).

% Copyright (c) 2018 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% This code was contributed by Stephen Sangwine and Eckhard Hitzer. See:
%
% Stephen J. Sangwine and Eckhard Hitzer,
% Polar decomposition of complexified quaternions and octonions, School of
% Computer Science and Electronic Engineering, University of Essex, UK
% Technical Report CES-535, April 2019. ISSN: 1744-8050 (to be submitted to
% Advances in Applied Clifford Algebras).

narginchk(1, 2), nargoutchk(2, 3)

if nargin == 1
    L = 'R'; % Supply the default value.
else
    if ~isa(L, 'char')
        error('Second parameter must be a character (L or R)')
    end
end

if L ~= 'L' && L ~= 'R'
    error('L must have the value ''L'' or ''R''.');
end

if nargout == 3
    % We must compute n, the semi-norm, and a normalised biquaternion to
    % work on from here on.
    
    n = abs(q);
    b = unit(q);
else
    if any(abs(abs(q(:))) - 1 > 100 * eps)
        warning(['Input parameter does not have unit modulus, third', ...
                 ' output parameter required.'])
    end
    b = q;
end

t = unit(real(b));

if strcmp(L, 'L')
    h = b .* conj(t);
else
    h = conj(t) .* b;
end

end

% $Id: polar.m 1025 2019-04-18 14:22:37Z sangwine $