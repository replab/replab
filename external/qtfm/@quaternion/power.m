function Z = power(X, Y)
% .^   Array power.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2006, 2008, 2016 Stephen J. Sangwine
%                                      and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

% This function can handle left and right parameters of the same size, or
% cases where one or other is a scalar.  The general case of a matrix or
% vector raised to a matrix or vector power requires elementwise operations
% and is handled using a general formula using logarithms, even though some
% of the elements of the right argument may be special cases (discussed
% below).

% When the right operand is a scalar, some special cases are handled using
% specific formulae because of the greater accuracy or better speed
% available. E.g. for Y == -1, the elementwise inverse is used, for Y == 2,
% elementwise squaring is used.

% For a power of ± 1/2, the sqrt function is used, with or without a
% reciprocal.

if isscalar(Y)
    
    % Y is a scalar. Check for and handle the various powers that are dealt
    % with as special cases.
    
    if     Y == -2, Z = (X .* X) .^ -1; % Use the next case recursively.
    elseif Y == -1
        Z = conj(X) ./ normq(X); % I.e. elementwise inverse. If X has
                                 % zero norm this will give a NaN.
    elseif Y == 0
        Z = ones(size(X));
    elseif Y == 1
        Z = X;
    elseif Y == 2
        Z = X .* X;
    elseif Y == 1/2
        Z = sqrt(X);
    elseif Y == -1/2
        Z = sqrt(X .^ -1); % Use the case Y == -1 above recursively.
    else
        Z = general_case(X, Y);
    end
    
elseif isscalar(X)

    % X is a scalar, but Y is not (otherwise it would have been handled
    % above). The general case code will handle this, since it will expand
    % X to the same size as Y before pointwise multiplication.
    
    Z = general_case(X, Y);

else
    
    % Neither X nor Y is a scalar, therefore we have to use the general
    % method. This will work only if the sizes are compatible. From Matlab
    % R2016b, 'compatible' has a looser interpretation based on implicit
    % singleton expansion. Rather than try to check for compatibility here,
    % we have removed the size check and we leave it to the Matlab code to
    % raise and error if the sizes are incompatible.
        
    Z = general_case(X, Y);

end

function Z = general_case(X, Y)
% The formula used here is taken from:
%
% A quaternion algebra tool set, Doug Sweetser,
%
% http://www.theworld.com/~sweetser/quaternions/intro/tools/tools.html
        
Z = exp(log(X) .* Y); % NB log(X) is the natural logarithm of X.
                      % (Matlab convention.)

% $Id: power.m 1004 2017-11-15 17:14:09Z sangwine $
