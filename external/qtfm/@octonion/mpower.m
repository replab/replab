function Z = mpower(X, Y)
% ^   Matrix power.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Based on the corresponding quaternion function.

narginchk(2, 2), nargoutchk(0, 1)

% There are three cases:
%
% 1. X and Y are both scalar. Handled by the .^ (power) function.
% 2. X is a square quaternion matrix and Y is a scalar. For Y == 2 it can
%    be computed by repeated multiplication. Otherwise non-associativity
%    makes the result ambiguous and we raise an error.
% 3. Neither of X or Y is scalar. This is an error, just as in the Matlab ^
%    function. There is no way to raise an (octonion) matrix to an
%    (octonion) matrix power.

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
        elseif Y > 2;
            error(['Powers of an octonion matrix are ambiguous due to ' ...
                   'non-associativity, not implemented'])
            % A^3 could mean A * (A^2) or (A^2) * A and the results would
            % differ.
        end
    end
end

help quaternion/mpower;
error('Cannot handle the given parameters yet.')

% $Id: mpower.m 1004 2017-11-15 17:14:09Z sangwine $
