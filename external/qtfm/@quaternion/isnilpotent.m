function tf = isnilpotent(q, tol)
% ISNILPOTENT  True where any element of q is a nilpotent to within the
% tolerance given (optionally) by the second parameter.

% Copyright (c) 2019 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 1
    tol = 4 .* eps; % The tolerance was not specified, supply a default.
end

% Reference:
%
% Stephen J. Sangwine and Daniel Alfsmann,
% Determination of the Biquaternion Divisors of Zero, Including the
% Idempotents and Nilpotents, Advances in Applied Clifford Algebras, 20,
% (2010), 401–410. DOI 10.1007/s00006-010-0202-3.

% Theorem 5 of the above paper gives the conditions for a biquaternion to
% be a nilpotent. The biquaternion must be pure, and the norm must be zero.

% We also need to check for a non-zero imaginary part, otherwise we would
% return true for a value of zero, which is not a nilpotent in the sense
% needed here.

if ispure(q)
    tf = abs(normq(imag(q))) > tol & ...
         abs(abs(normq(q))) < tol;
else
    % Nilpotents must be pure, so we check that the scalar part is close to
    % zero, then check the norm of the vector part is also close to zero.
    
    tf = abs(normq(imag(q))) > tol & ...
         s(q) < tol & abs(abs(normq(v(q)))) < tol;
end

end

% $Id: isnilpotent.m 1031 2019-04-20 12:49:25Z sangwine $
