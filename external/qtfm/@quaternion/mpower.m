function Z = mpower(X, Y)
% ^   Matrix power.
% (Quaternion overloading of standard Matlab function.)
% For cases that can be computed using .^ (power), this function calls the
% power function. Other cases are not yet implemented for quaternions but
% may be added at a later date. Please contact the authors of QTFM.

% Copyright (c) 2008, 2009 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

% There are three cases:
%
% 1. X and Y are both scalar. Handled by the .^ (power) function.
% 2. X is a square quaternion matrix and Y is a scalar. For small integer
%    powers this can be computed by repeated multiplication. Otherwise it
%    is not yet implemented.
% 3. Neither of X or Y is scalar. This is an error, just as in the Matlab ^
%    function. There is no way to raise a (quaternion) matrix to a
%    (quaternion) matrix power.

if isscalar(X) && isscalar(Y) % Case 1.
    Z = power(X, Y); % The power function will check the sanity of X and Y.
    return;
end

if ~isscalar(X) && ~isscalar(Y) % Case 3.
   error('At least one operand must be scalar.') 
end

if isscalar(Y) && isnumeric(Y) % Possibly case 2,
    if length(size(X)) == 2 && size(X, 1) == size(X, 2) % if X is square.
        % Case 2. We can handle small positive integer values by repeated
        % multiplication. We do this only for certain values because errors
        % will build up.
        if Y == 0
            Z = eyeq(size(X)); % Matlab handles this one, so why not?
            return;
        elseif Y == 1
            Z = X; return;
        elseif Y == 2;
            Z = X * X; return;
        elseif Y == 3;
            Z = X * X * X; return;
        % We choose not to go beyond 3. Would it be better to compute
        % X * X * X * X for Y == 4, or would it be better to compute X * X
        % and then square it? What about Y == 6? Where to stop? Decided to
        % leave this for the more general case.
        end
    end
end

help quaternion/mpower;
error('Cannot handle the given parameters yet.')

% TODO Implement the function mpower for the general case. Does the
% Cayley-Hamilton Theorem apply to quaternion matrices? This could provide
% a general and efficient way to compute higher matrix powers.

% $Id: mpower.m 1024 2019-04-18 14:20:12Z sangwine $
