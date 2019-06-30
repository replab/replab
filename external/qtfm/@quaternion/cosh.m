function Y = cosh(X)
% COSH   Hyperbolic cosine.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(X)
    
    % X is a real quaternion, and we compute the hyperbolic cosine of an
    % isomorphic complex number using the standard Matlab cosh function,
    % then construct a quaternion with the same axis as the original
    % quaternion.
    
    Y = isoquaternion(cosh(isocomplex(X)), X);
else
    
    % X is a complex quaternion, and therefore we cannot use the method
    % above for real quaternions, because it is not possible to construct
    % an isomorphic complex number. We use instead a fundamental formula
    % for the hyperbolic cosine.
    
    Y = (exp(X) + exp(-X)) ./ 2;
end

% $Id: cosh.m 1008 2018-05-09 20:52:12Z sangwine $
