function theta = angle(q, a)
% ANGLE or phase of octonion.
% (Octonion overloading of standard Matlab function.)
%
% The second parameter is optional. If omitted, the result will always be
% in the range (0, pi), that is, the octonion q is regarded as being in
% the upper half of a complex plane defined by the axis of q. A reference
% unit vector defining the direction of the positive imaginary axis of q
% may be supplied as the second parameter, in which case the result may
% be in the full range from -pi to +pi. This reference unit vector is used
% to define the north pole of a hemisphere of 7-space, so that if the axis
% of q lies in this hemisphere, the angle is in the range (0,pi). If the
% axis of q lies in the other (southern) hemisphere, then the angle
% returned will lie in (pi, 2*pi).
%
% The optional argument must be either the same size as q or be a scalar
% (in which case the same value is used for all elements of q).

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 2
    if ~isa(a, 'octonion')
        error('Second argument must be an octonion (array).');
    end
    if ~isreal(a)
        error('Second argument must be a real octonion (array).');
    end
    if ~isempty(a.w) || any(any(abs(abs(a) - 1) > eps))
        error('Second argument must be a pure unit octonion (array).');
    end
    if ~(all(size(q) == size(a)) || isscalar(a))
        error(['The two parameters must be the same size,' ...
               ' or the second one must be a scalar.']);
    end
end

if ispure(q)
    
    % q is a pure octonion, with no scalar part, and by definition, the
    % angle is a right-angle, even for complexified octonions. Return a
    % matrix of the same size as q, with pi/2 at every position.
    % Rationale : we can construct an exact result in this case with
    % minimal calculation, whereas if we use the algorithm below for full
    % octonions, the result may be subject to small errors.
    
    theta = ones(size(q), class(part(q, 2))) .* (pi ./ 2);
    
else
    
    % q is a full octonion. We adopt two methods here, depending on
    % whether q is a real octonion or a complexified octonion.
 
    if isreal(q)

        % All elements of q are real. Therefore we can use the simple
        % method of constructing a complex value from q and using the
        % Matlab angle function to compute the angle. We assume that the
        % Matlab angle function does a good job on complex numbers with
        % small modulus, and do not try to better it.

        if nargin == 2
            theta = angle(isocomplex(q, a));
        else
            theta = angle(isocomplex(q));
        end
    else
        
        % One or more elements of q is/are complex. We cannot construct an
        % isomorphic complex array in this case, so we have to resort to a
        % more fundamental method using the complex scalar and vector parts
        % of q. See the appendix below for details of the derivation of the
        % formula used. TODO This code needs to be improved to handle cases
        % with small modulus properly and return a zero angle.

        x  =       s(q);
        y2 = normo(v(q)); % This is equivalent to y.^2, y = abs(v(q))
       
        theta = - 1i .* log((x + 1i .* sqrt(y2)) ./ sqrt(x.^2 + y2));
    end
end

% Appendix: calculation of arctangent(y, x) when y and x are complex.
%
% (See the quaternion version of this function for details.)

% $Id: angle.m 1004 2017-11-15 17:14:09Z sangwine $
