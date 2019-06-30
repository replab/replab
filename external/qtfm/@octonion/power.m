function Z = power(X, Y)
% .^   Array power.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
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
    
    % Y is a scalar.
    
    if Y == fix(Y)

        % Y is an integer, we can handle some cases.
        
        if Y == -2
            Z = (X .* X) .^ -1; % Use the next case recursively.
        elseif Y == -1
            Z = conj(X) ./ normo(X); % I.e. elementwise inverse. If X has
            % zero norm this will give a NaN.
        elseif Y == 0
            Z = ones(size(X));
        elseif Y == 1
            Z = X;
        elseif Y == 2
            Z = X .* X;
        else
            % None of the above fits, but for small integer powers we can
            % use repeated multiplication (the octonions are power
            % associative).
            if abs(Y) < 10 % This limit is arbitrary. 
                Z = X;
                for j = 1:abs(Y)-1 % NB We have excluded |Y|< 3 already.
                    Z = Z .* X;
                end
                if Y < 0
                   % The power is negative, we need to take the inverse.
                   Z = Z .^ -1;
                end
            else
                Z = general_case(X, Y); % Give up and resort to logs.
            end
        end
    else
        
        % Y is not an integer, but there may be some special cases.
    
        if Y == 1/2
            Z = sqrt(X);
        elseif Y == -1/2
            Z = sqrt(X .^ -1); % Use the case Y == -1 above recursively.
        elseif Y ~= fix(Y)
            % Y is not an integer, therefore the general case has to apply.
            Z = general_case(X, Y);
        else
            
        end
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

% Whether it works for octonions has not been checked, hence the warning!
warning('Octonion elementwise power has not yet been checked.');
        
Z = exp(log(X) .* Y); % NB log(X) is the natural logarithm of X.
                      % (Matlab convention.)

% $Id: power.m 1004 2017-11-15 17:14:09Z sangwine $
