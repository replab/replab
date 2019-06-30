function [p, m] = opd(q, a, b)
% OPD Orthogonal plane decomposition of a quaternion.
%
% This function decomposes a quaternion (array) q element-by-element into
% components defined by the scalar quaternions a and b. The third parameter
% is optional: if omitted it is set equal to the second. The geometric
% effect of the decomposition is as follows:
%
% 1. If a is a pure quaternion, and b is omitted, p is parallel to a, and 
%    m is in the plane normal to a. (Parallel/perpendicular decomposition.) 
% 2. If a is a full quaternion, and b is omitted, ....
% 3. If a and b are pure quaternions, ....
% 4. If a and b are full quaternions, ....
% Other possibilities may exist.

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 3), nargoutchk(2, 2) 

% References:
%
% E. Hitzer and S. J. Sangwine,
% 'The orthogonal planes split of quaternions', in: K. Guerlebeck (ed.),
% 9th International Conference on Clifford Algebras and their Applications
% in Mathematical Physics, Weimar, Germany, 15-20 July 2011.
%
% Ell, T. A. and Sangwine, S. J.,
% 'Quaternion Involutions and Anti-Involutions',
% Computers and Mathematics with Applications, 53(1), January 2007, 137-143.
% doi:10.1016/j.camwa.2006.10.029.
%
% H. S. M. Coxeter, 'Quaternions and reflections',
% American Mathematical Monthly, 53(3), 136?146, 1946.

if ~isa(q, 'quaternion')
    error('First parameter must be a quaternion.')
end

if ~isa(a, 'quaternion')
    error('Second parameter must be a quaternion.')
end

if ~isscalar(a)
    error('Second parameter must be a scalar quaternion.')
end

if nargin == 3
    if ~isa(b, 'quaternion')
        error('Third parameter must be a quaternion.')
    end
    
    if ~isscalar(a)
        error('Second parameter must be a scalar quaternion.')
    end
else
    b = a;
end

% The scalar parameters a and b are required to have unit modulus for the
% geometric operations discussed above to be correct. Rather than ask this
% of the user and check for unit modulus, we accept any modulus, but
% unitize the quaternions here. If they already have unit modulus, this
% operation has no ill effect.

t = unit(a) .* q .* unit(b);

p = (q + t)./2; % p stands for plus.
m = (q - t)./2; % m stands for minus.

% TODO If a is a pure quaternion, and b is omitted, the above formulae
% result by definition in a pure quaternion result for p, and we should
% remove the scalar part before return.

end

% $Id: opd.m 1004 2017-11-15 17:14:09Z sangwine $
