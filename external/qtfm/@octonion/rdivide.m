function r = rdivide(a, b)
% ./  Right array divide.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

% There are three cases to be handled:
%
% 1. Left and right arguments are octonions.
% 2. The left argument is an octonion, the right is not.
% 3. The right argument is an octonion, the left is not.
%
% In fact, cases 1 and 3 can be handled by the same code. Case 2
% requires different handling.

if isa(b, 'octonion')
    
    % The right argument is an octonion. We can handle this case by
    % forming its elementwise inverse and then multiplying. Of course,
    % if any elements have zero norm, this will result in NaNs.
     
    r = a .* b .^ -1;
    
else
    
    % The right argument is not an octonion. We assume therefore
    % that if we divide components of the left argument by the right
    % argument, that Matlab will do the rest. Obviously if the right
    % argument is zero, there will be a divide by zero error.
    % TODO This is somewhat dangerous - what if b is a quaternion?

    r = overload(mfilename, a, b);
end

% $Id: rdivide.m 1004 2017-11-15 17:14:09Z sangwine $
