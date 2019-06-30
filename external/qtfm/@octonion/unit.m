function r = unit(a)
% UNIT octonion. Divides an octonion by its own modulus.
% The result is an octonion with unit modulus.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

m = abs(a);

% Since m could be complex, the warning check below has to
% use abs again to get a real result for comparison with eps.

if any(any(abs(m) < eps))
    warning('QTFM:information', ...
            ['At least one element has zero modulus, '...
             'and divide by zero will occur.']);
end

r = a ./ m;
 
% Dividing by a small modulus can result in numerical errors such that the
% result does not have unit modulus. This is especially likely with complex
% octonions, so we perform a check here.

% The modulus (and norm) of each element of r should be either 1 or 1 + 0i. 
% In either case, subtracting ones should result in a real or complex
% number with very small real and imaginary parts, which we compare with
% epsilon using an arbitrary scale factor which we choose so that the test
% is not too sensitive.

d = normo(r) - 1;
n = nnz(abs(real(d)) > 1e3.*eps | abs(imag(d)) > 1e3.*eps);

if n > 0
    if n == 1
        warning('QTFM:information', ...
            ['One element of the result of the unit '...
            'function has a modulus which is not accurately 1.']);
    else
        warning('QTFM:information', ...
            [num2str(n), ' elements of the result of the unit '...
            'function have a modulus which is not accurately 1.']);
    end
end

% Discussion note. See the note in the corresponding quaternion function
% for possible issues with accuracy of unit complex octonions.

% $Id: unit.m 1004 2017-11-15 17:14:09Z sangwine $

