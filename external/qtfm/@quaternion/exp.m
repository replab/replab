function Y = exp(X)
% EXP  Exponential
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2019 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Modified 21 September 2005 to handle undefined axis correctly.

narginchk(1, 1), nargoutchk(0, 1)

if isempty(X.w)
    
    % X is a pure quaternion (array).  We are computing exp(X) for each
    % element of the array, and we need to factor X to give a unit pure
    % quaternion, a, and a real or complex 'angle' denoted below by theta.
    % In fact this angle is simply the absolute value of X. We then have
    % X = a .* theta, and we can calculate exp(X) from cos(theta) +
    % a .* sin(theta) (de Moivre's formula).
    %
    % If the modulus of any element of X is zero or very small, the axis
    % (a) will be undefined, and we have to take care to handle this
    % without returning NaNs, since although the axis is undefined, the
    % exponential is not. A second reason why the modulus may be small is
    % that the element in question is a nilpotent. We test for this
    % separately, because in this case the exponential is 1 + X, as can be
    % seen from the power series for e^x.
    
    theta     = abs(X);           % Since theta could be complex, we need
    undefined = abs(theta) < eps; % abs() again here to compare with eps.
    
    cc = cos(theta); % If theta is complex, the cosine and sine will yield
    ss = sin(theta); % cosh and I.*sinh of theta.

    theta(undefined) = 1; % This prevents divide by zero in the element
                          % positions where the axis is undefined.

    a = X ./ theta; % This is equivalent to axis(X) but without warnings
                    % about undefined axis because values of theta less
                    % than eps have been replaced with 1. This means the
                    % values at positions corresponding to undefined axis
                    % are incorrect, since we have copied the (small) value
                    % from X unmodified. However, these values will be
                    % multiplied by a zero or very small (less than eps)
                    % sine value below, so are harmless.

    Y = cc + a .* ss;
    
    % Now deal with nilpotent elements if any.
    
    if any(undefined(:))
        % Note on subscripted references: these do not work inside a class
        % method (which this is). See the file 'implementation notes.txt',
        % item 8.
        
        N = false(size(X));

        % N(undefined) = isnilpotent(X(undefined))
        
        N(undefined) = isnilpotent( ...                           
                       subsref(X, substruct('()', {undefined})));
        
        if any(N(:)) % There may be no nilpotents, in which case, skip.
            Y = subsasgn(Y, substruct('()', {N}), ...  % Y(N) = ...
                1 + subsref(X, substruct('()', {N}))); %        1 + X(N);
        end
    end
else
    
    % X is a full quaternion (array).  We use a recursive call to compute
    % the exponential of the vector part, and the standard Matlab function
    % to compute the exponential of the scalar part (which may be complex).
    % The result is the elementwise product of the two.

    Y = exp(s(X)) .* exp(v(X));
end

% $Id: exp.m 1031 2019-04-20 12:49:25Z sangwine $
