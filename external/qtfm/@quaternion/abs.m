function a = abs(q)
% ABS Absolute value, or modulus, of a quaternion.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(q)
    
    % We use here a method based on Cayley-Dickson form and the Matlab
    % function hypot, which avoids overflow for large values.
    
    if isempty(q.w)
        a = hypot(abs(q.x), abs(complex(q.y, q.z)));
    else
        a = hypot(abs(complex(q.w, q.x)), abs(complex(q.y, q.z)));
    end
    
else
    
    % Here, for biquaternions, we cannot use hypot and we use the obvious
    % formula, which was used for all quaternions up to version 1.4 of this
    % file.
    
    if isempty(q.w)
        a = sqrt(         q.x.^2 + q.y.^2 + q.z.^2);
    else
        a = sqrt(q.w.^2 + q.x.^2 + q.y.^2 + q.z.^2);
    end
    
end

end

% $Id: abs.m 1004 2017-11-15 17:14:09Z sangwine $
