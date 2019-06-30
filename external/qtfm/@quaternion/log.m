function Y = log(X)
% LOG    Natural logarithm.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2006, 2010, 2019
%               Stephen J. Sangwine and Nicolas Le Bihan.
%
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(X)
    
    % X is a real quaternion, and we compute the logarithm of an isomorphic
    % complex number using the standard Matlab log function, then construct
    % a quaternion with the same axis as the original quaternion.
    
    Y = isoquaternion(log(isocomplex(X)), X);
else
    
    % X is a complex quaternion (or at least some elements of X are), and
    % therefore we cannot use the method above for real quaternions,
    % because it is not possible to construct an isomorphic complex number.
    % The algorithm used is documented in the appendix. However, this will
    % not work for the exponential of a nilpotent, and special treatment is
    % needed if these are present. See the exp function for background.
    
    N = (abs(abs(scalar(X) - 1)) < eps) & ...
        isnilpotent(v(X)); % This will detect values of the form
                           % 1 + X = exp(X) where X is a nilpotent.
                                       
    % TODO Replace the scalar parts of any values indexed by N with 0.
    % Do the calculation below only on the other elments of N. Or find a
    % better algorithm for the log that will work in all cases?

    if any(N(:))
        warning(['At least one array element is the exponential of a ', ...
                 'nilpotent, logarithm value will be incorrect.']);
    end
    
    Y = quaternion(log(normq(X))./2, axis(X) .* angle(X));
     
end

% Appendix.
%
% The calculation of log(X) is not difficult to derive, as follows.
%
% First note that exp(Y) = X (definition of logarithm). Then write
% X in polar form as r .* exp(mu .* theta). Then we have:
%
% log(X) = Y = log(r) + log(exp(mu .* theta)) = log(r) + mu .* theta.
% 
% Rather than calculate r = abs(X), we use normq because this
% avoids calculating a square root. We compensate by halving log(r).
% mu and theta are obtained using the axis and angle functions. Note
% that we don't need to add the two results, since the first is a
% real or complex value, and the second is a pure quaternion (real
% or complex). Therefore we can compose them using the quaternion
% constructor function.

% $Id: log.m 1032 2019-04-20 14:20:06Z sangwine $
