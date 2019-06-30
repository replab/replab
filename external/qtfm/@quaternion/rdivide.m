function r = rdivide(a, b)
% ./  Right array divide.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2006, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

% There are three cases to be handled:
%
% 1. Left and right arguments are quaternions.
% 2. The left argument is a quaternion, the right is not.
% 3. The right argument is a quaternion, the left is not.
%
% In fact, cases 1 and 3 can be handled by the same code. Case 2
% requires different handling.

if isa(b, 'quaternion')
    
    % The right argument is a quaternion. We can handle this case by
    % forming its elementwise inverse and then multiplying. Of course,
    % if any elements have zero norm, this will result in NaNs.
     
    r = a .* b .^ -1; % Changed 24/2/2006 to use .^ -1 instead of inv.
    
else
    
    % The right argument is not a quaternion. We assume therefore
    % that if we divide components of the left argument by the right
    % argument, that Matlab will do the rest. Obviously if the right
    % argument is zero, there will be a divide by zero error.
    
    r = a;
    if ~isempty(a.w)
        r.w = r.w ./ b;
    end
    r.x = a.x ./ b;
    r.y = a.y ./ b;
    r.z = a.z ./ b;
end

% $Id: rdivide.m 1004 2017-11-15 17:14:09Z sangwine $
