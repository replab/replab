function [T, N1, N2, N3, J, K, L, M] = frenet(q)
% FRENET   Frenet-Serret frames of a quaternion sequence or curve.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO Consider an alternative formulation as set out by Kurt Nalty (see
% references below) in which the curvature, torsion and bitorsion are
% vectors). This could be compared with the current implementation, and
% thus checked for correctness before a decision to release it is made.

% Given a vector of quaternions (full or pure), computes the vectors of
% tangents, normals, binormals, and trinormals and the curvatures. The
% output parameters are as follows. The first four results are unit
% quaternions; the remaining four are scalars.
%
% T  are the tangents;
% N1 are the normals;
% N2 are the binormals;
% N3 are the trinormals;
%
% J is the moduli of the tangents, so that non-unit tangents are J .* T;
% K is the curvature; again K .* N1 gives un-normalised normals;
% L is the second curvature or torsion;
% M is the third curvature or bitorsion.
%
% Note that the first and last values of T are not necessarily meaningful
% because T is an approximation to the derivative, and therefore cannot be
% correctly computed for the first and last values in q. The same applies
% to N1, N2 and N3, but as these are computed using successive derivatives,
% more values at the extremities are not necessarily meaningful.

% References:
%
% K. Bharathi and M. Nagaraj,
% 'Quaternion Valued Function of a Real Variable Serret-Frenet Formula',
% Indian J. Pure and Applied Mathematics, 18(6), 507-511, June 1987.
%
% 'Frenet-Serret formulas', subsection: 'Formulas in n dimensions',
% Wikipedia (English), accessed 27 April 2011.
% http://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas
%
% Kurt Nalty, 'Trajectories and curves in 3 and 4 dimensions',
% self-published article, 11pp, 19 January 2007, available at
% http://kurtnalty.com/.
%
% Comments: the article by Bharathi and Nagaraj was probably the first to
% explain the Frenet-Serret frames in 4-space using quaternions, but the
% Wikipedia article makes clear that Camille Jordan worked out the
% n-dimensional case in 1874. The article by Kurt Nalty is of value for its
% clear tutorial explanation of the formulae.

narginchk(1, 1), nargoutchk(1, 8)

if ~isvector(q)
    error('The input parameter must be a quaternion vector.')
end

% TODO are there other error checks that should be done?

[J, T ] = absunit(derivative(q));
[K, N1] = absunit(derivative(T));
[L, N2] = absunit(derivative(N1) + K .* T);
[t, N3] = absunit(derivative(N2) + L .* N1); M = t - K;

end

function [a, u] = absunit(x)
% Splits a quaternion into its modulus and a unit quaternion.
a = abs(x); u = unit(x);
end

% TODO Consider whether to allow a 9th parameter to control the derivative
% algorithm. We could for example use diff to compute the derivative, but
% this has the disadvantage that the result gets shorter by one element at
% each derivative step.

function y = derivative(x)
% Compute the derivative of x. This is done below using a central finite
% difference. Since it uses circular shift there is no change in the length
% of the result, but this does mean that the first and last values are
% meaningful only if the sequence in x is periodic and contains an integral
% number of periods (informally, x is cyclic continuous).

y = circshift(x, [1,-1]) - circshift(x, [1,+1]);
end

% $Id: frenet.m 1004 2017-11-15 17:14:09Z sangwine $
