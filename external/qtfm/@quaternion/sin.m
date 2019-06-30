function Y = sin(X)
% SIN    Sine.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(X)
    
    % X is a real quaternion, and we compute the sine of an isomorphic
    % complex number using the standard Matlab sin function, then
    % construct a quaternion with the same axis as the original quaternion.
    
    Y = isoquaternion(sin(isocomplex(X)), X);
else
    
    % X is a complex quaternion, and therefore we cannot use the method
    % above for real quaternions, because it is not possible to construct
    % an isomorphic complex number. We use instead a fundamental formula
    % for the sine in terms of the exponential.
    
    Y = (exp(1i .* X) - exp(-1i .* X)) ./ (2 .* 1i);
end;

% $Id: sin.m 1004 2017-11-15 17:14:09Z sangwine $
