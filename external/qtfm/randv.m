function S = randv(varargin)
% RANDV   Creates uniformly distributed unit vectors (pure quaternions).
% Accepts the same parameters as the Matlab function rand (q.v.). The
% result may be scalar, vector, matrix or array depending on the parameters
% supplied. Each pure quaternion returned is the result of two calls on
% rand, and hence randv modifies the state of the generator used by rand.
% To initialise the generator or control the choice of generator, use rand.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan and Salem Said.
% See the file : Copyright.m for further details.

nargoutchk(0, 1)

% Reference:
%
% James Arvo, 'Random rotation matrices', Graphics Gems II, James Arvo (ed)
% Academic Press, 1991, pp355-6.
%
% [The algorithm used here is only part of Arvo's algorithm, to generate a
% random axis in 3-space for a rotation matrix. It is also a simplification
% of the algorithm by Shoemake used in randq.m (q.v.), which is clearly
% inspired by Arvo's algorithm (and corrects an error in it).]

z = rand(varargin{:}) .* 2 - 1; % z is uniformly distributed in (-1, +1).
r = sqrt(1 - z.^2);

theta = rand(varargin{:}) .* 2 .* pi; % Uniformly distributed in (0, 2*pi).

S = quaternion(r .* cos(theta), r .* sin(theta), z);

end

% Note: a simpler, but slower algorithm (by about a factor of 1.8) is the
% following, using three calls on randn instead of two on rand. Because the
% three components of the quaternion are normally distributed, the
% directions are uniformly distributed in 3-space. This seems to be fairly
% well-known, but it can be verified in the following reference, section 5,
% p209:
%
% J. K. Mackenzie and M. J. Thomson, "Some statistics associated with the
% random disorientation of cubes", Biometrika, 44(1/2), June 1957, 205-210.
%
% S = unit(quaternion(randn(varargin{:}), ...
%                     randn(varargin{:}), ...
%                     randn(varargin{:})));

% $Id: randv.m 1004 2017-11-15 17:14:09Z sangwine $
