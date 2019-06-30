function S = randq(varargin)
% RANDQ   Creates uniformly distributed unit quaternions (points on the
% 4-sphere), or uniformly distributed rotations in 3-space represented in
% their quaternion form. Accepts the same parameters as the Matlab function
% rand (q.v.). The result may be scalar, vector, matrix or array depending
% on the parameters supplied. Each quaternion returned is the result of
% three calls on rand, and hence randq modifies the state of the generator
% used by rand. To initialise the generator or control the choice of
% generator, use rand.
%
% The result will have uniformly distributed axes and Gaussian-distributed
% angles. This rather non-intuitive result is discussed by Shoemake (see
% reference). If the quaternions supplied by this function are applied to a
% constant pure quaternion (e.g. qi) using the formula conj(S) .* qi .* S,
% the result will be randomly oriented pure quaternions with a uniform
% distribution in 3-space: this is what is meant by 'uniformly distributed
% rotations (but note that this result can be obtained more directly with
% randv (q.v.)).

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan and Salem Said.
% See the file : Copyright.m for further details.

% Reference:
%
% Ken Shoemake, 'Uniform random rotations', in David Kirk (ed), Graphics
% Gems III, Academic Press, 1992, pp124-132.
%
% [The algorithm used here is described on pages 129-130.]

nargoutchk(0, 1)

x0 = rand(varargin{:});

theta1 = rand(varargin{:}) .* 2 .* pi;
theta2 = rand(varargin{:}) .* 2 .* pi;

s1 = sin(theta1); s2 = sin(theta2);
c1 = cos(theta1); c2 = cos(theta2);

r1 = sqrt(1 - x0);
r2 = sqrt(x0);

S = quaternion(s1 .* r1, c1 .* r1, s2 .* r2, c2 .* r2);

end

% Note: a more concise, but slower (by about 20%) algorithm is given below:
% 
% S = unit(quaternion(randn(varargin{:}), ...
%                     randn(varargin{:}), ...
%                     randn(varargin{:}), ...
%                     randn(varargin{:})));
%
% This algorithm is discussed in Shoemake's article (above), on page 129.
% That it is not obvious is shown by the observation that the corresponding
% algorithm in 3-dimensions gives a different distribution of the
% components after normalisation (uniform), whereas here the components are
% neither uniform, nor Gaussian-distributed after normalisation. However,
% what can be verified easily is that normalising a quaternion with
% Gaussian distributed components gives the same component distributions
% as the Shoemake algorithm used above.

% $Id: randq.m 1004 2017-11-15 17:14:09Z sangwine $
