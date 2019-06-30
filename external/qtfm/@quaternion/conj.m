function c = conj(a, F)
% CONJ   Quaternion conjugate.
% (Quaternion overloading of Matlab standard function.)
%
% Implements three different conjugates:
%
% conj(X) or
% conj(X, 'hamilton') returns the quaternion conjugate.
% conj(X, 'complex')  returns the complex conjugate.
% conj(X, 'total')    returns the 'total' conjugate equivalent to
%                     conj(conj(X, 'complex'), 'hamilton')

% Copyright (c) 2005, 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

c = a; % Copy the input parameter to make the return result.

if nargin == 1
    % The most common case, just compute the conjugate and we are done.
    c.x = -c.x; c.y = -c.y; c.z = -c.z;
else
    if ~strcmp(F, 'hamilton') && ~strcmp(F, 'complex') && ~strcmp(F, 'total')
        error('Second parameter value must be ''hamilton'', ''complex'' or ''total''.')
    end
    
    switch F
        case 'hamilton'
            c.x = -c.x; c.y = -c.y; c.z = -c.z;
        case 'complex'
            c = overload(mfilename, a); % Conjugate the components of a.
        case 'total'
            c = conj(overload(mfilename, a));
        otherwise
            error('Bad value for second parameter.');
    end
end

end

% $Id: conj.m 1004 2017-11-15 17:14:09Z sangwine $
