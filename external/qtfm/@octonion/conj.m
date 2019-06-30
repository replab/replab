function c = conj(a, F)
% CONJ   Octonion conjugate.
% (Octonion overloading of Matlab standard function.)
%
% Implements three different conjugates:
%
% conj(X) or
% conj(X, 'cayley')  returns the octonion conjugate.
% conj(X, 'complex') returns the complex conjugate.
% conj(X, 'total')   returns the 'total' conjugate equivalent to
%                    conj(conj(X, 'complex'), 'cayley')

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

c = a; % Copy the parameter in order to create an octonion result.
    
if nargin == 1
    % Deal with the most common case and return.
    c.a = conj(c.a); % The quaternion conjugate, of course, which
                     % negates all but the scalar part.
    c.b = - c.b;     % Negate the other four components.
else

    if ~strcmp(F, 'cayley') && ~strcmp(F, 'complex') && ~strcmp(F, 'total')
        error('Second parameter value must be ''cayley'', ''complex'' or ''total''.')
    end
    
    switch F
        case 'cayley'
            c.a = conj(c.a); % The quaternion conjugate, of course, which
                             % negates all but the scalar part.
            c.b = - c.b;     % Negate the other four components.
        case 'complex'
            c.a = conj(c.a, F); % Take the complex conjugate of each quaternion
            c.b = conj(c.b, F); % component.
        case 'total'
            c = conj(conj(a, 'complex'), 'cayley');
        otherwise
            error('Bad value for second parameter.');
    end
end

end

% $Id: conj.m 1004 2017-11-15 17:14:09Z sangwine $
