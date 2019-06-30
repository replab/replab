function r = le(a, b)
% <=  Less than or equal
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% The Matlab function operating on complex values compares only the real
% parts, so here we do a comparison of the scalar parts. If these are
% complex the Matlab function will ignore the imaginary parts.

narginchk(2, 2), nargoutchk(0, 1)

if isa(a, 'quaternion') && isa(b, 'quaternion')
    r = scalar(a) <= scalar(b);
else
    % One of the arguments is not a quaternion (the other must be, or
    % Matlab would not call this function). The non-quaternion argument
    % must be a numeric (if we don't impose this restriction it could be
    % anything such as a cell array or string which makes no sense at all
    % to compare with a quaternion).
            
    if isa(a, 'quaternion') && isa(b, 'numeric')
        r = scalar(a) <= b;
    elseif isa(a, 'numeric') && isa(b, 'quaternion')
        r = a <= scalar(b);
    else
        error('Cannot compare a quaternion with a non-numeric');    
    end
end

% $Id: le.m 1004 2017-11-15 17:14:09Z sangwine $
