function tf = isdivisor(q, tol)
% ISDIVISOR  True where any element of q is a divisor of zero to within the
% tolerance given (optionally) by the second parameter.

% Copyright (c) 2019 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 1
    tol = eps; % The tolerance was not specified, supply a default.
end

% Reference:
%
% Stephen J. Sangwine and Daniel Alfsmann,
% Determination of the Biquaternion Divisors of Zero, Including the
% Idempotents and Nilpotents, Advances in Applied Clifford Algebras, 20,
% (2010), 401–410. DOI 10.1007/s00006-010-0202-3.

% Theorem 1 of the above paper gives the conditions for a quaternion to be
% a divisor of zero.

tf = abs(normq(real(q)) - normq(imag(q)))  < tol & ...
     abs(scalar_product(real(q), imag(q))) < tol;

end

% $Id: isdivisor.m 1031 2019-04-20 12:49:25Z sangwine $
