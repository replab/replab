function a = abs(o)
% ABS Absolute value, or modulus, of an octonion.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

if isreal(o)
    
    % We use here a method based on Cayley-Dickson form and the Matlab
    % function hypot, which avoids overflow for large values.
    
    a = hypot(abs(o.a), abs(o.b));    
else
    
    % Here, for bioctonions, we cannot use hypot.
    
    a = sqrt(normq(o.a) + normq(o.b)); 
end

end

% $Id: abs.m 1004 2017-11-15 17:14:09Z sangwine $
