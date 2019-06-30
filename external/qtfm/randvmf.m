function S = randvmf(mu, kappa, varargin)
% RANDVMF   Returns unit quaternions distributed on the 4-sphere according
% to a von Mises-Fisher distribution. See also RANDF for the 3-sphere case.
%
% The first parameter, mu, must be a full quaternion, defining the
% direction of the centre of the distribution in 4-space. The second
% parameter is the concentration parameter which controls the spread of the
% distribution on the sphere. It must be non-negative. A value of zero
% results in a uniform distribution on the sphere. Larger values result in
% greater concentration of the distribution in the mean direction mu.
%
% The remaining parameters are as for the Matlab function rand (q.v.).  The
% result may be scalar, vector, matrix or array depending on the parameters
% supplied.  Each pure quaternion returned is the result of at least two
% calls on rand, and two calls on randn, and hence randf modifies the state
% of the generator used by both rand and randn. To initialise the generator
% or control the choice of generator, use rand and/or randn.

% This code was rearranged in January 2016 (CVS revision 1.5) to vectorise
% the calculation of S. A side effect of this change is that calls to rand
% and randn are done in a different order to previous versions, which means
% that the function will not produce the same random values as it did, even
% if the random number generator has the same initial value.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% References:
%
% [The paper by Fisher is the definitive paper on the distribution on a
% sphere. Mardia and Jupp discuss it and the generalization to
% n-dimensions, but refer to Wood for details of how to generate a set of
% samples from the distribution. Wood's work is based on the paper by
% Ulrich. The report by Dhillon and Sra was used to code the implementation
% below, but they gave an algorithm for the general m-dimensional case. We
% have simplified the code for the case of three dimensions and also made a
% correction -- their Figure 4 does not make use of c after computing it.
% Wood's paper reveals that it should be used in the test on Z and U.]
%
% R. A. Fisher, "Dispersion on a sphere", Proceedings of the Royal Society
% of London Series A., 217, pp295-305, (1953).
%
% K. V. Mardia and P. E. Jupp, "Directional Statistics" (2nd edition), John
% Wiley (2000). ISBN 0-471-95333-4. [§9.3.]
%
% Gary Ulrich, "Computer Generation of Distributions on the m-Sphere",
% Applied Statistics, 33(2), pp158-163, (1984).
%
% A. T. A. Wood, "Simulation of the von-Mises Distribution", Communications
% in statistics : simulation & computation, 23, pp157-164, 1994.
%
% Inderjit S. Dhillon and Suvrit Sra, "Modeling Data using Directional
% Distributions", Technical Report TR-03-06, Department of Computer
% Sciences, The University of Texas at Austin, Austin, TX 78712, USA.
% 25 January, 2003.
% Accessed at: http://www.cs.utexas.edu/research/publications/ in May 2008.

nargoutchk(0, 1)

if ~isscalar(mu)
    error('The first parameter (mean direction) must be scalar.')
end

if ~isa(mu, 'quaternion') || ispure(mu)
    error('The first parameter (mean direction) must be a full quaternion.')
end

if ~isnumeric(kappa) || ~isscalar(kappa) || ~isreal(kappa)
    error('The second (concentration) parameter must be a numeric scalar.')
end

if kappa < 0
    error('The second (concentration) parameter must be non-negative.')
end

b = -kappa + sqrt(kappa.^2 + 1);

x0 = (1 - b) ./ (1 + b);

c = kappa .* x0 + 2 .* log(1 - x0.^2);

w = zeros(varargin{:}); % Initialise an array of the required size and
                        % shape, but we will address it in the loop using
                        % linear indexing for simplicity.

% TODO See comment in randf.m about vectorising the computation of w.

for i = 1:numel(w)
    while true
        z = beta; % Note 1.
        u = rand;
        t = (1 - (1 + b) .* z) ./ (1 - (1 - b) .* z);
        if kappa * t + 2 * log(1 - x0 .* t) - c >= log(u), break, end
    end
    w(i) = t;
end

% Now compute S in one vectorised step, using the values of w.

% We now need to create a point on the unit 3-sphere, uniformly
% distributed. This is done with the randv function. We then combine
% this value with w, keeping the modulus to unity, and we have the
% result we need, apart from a final rotation.

S = quaternion(w, sqrt(1 - w.^2) .* randv(size(w)));

% S now has elements with the vMF distribution, but mean direction
% (1,0,0,0) i.e. aligned with the scalar axis. We must now rotate the
% elements of S by the rotation that takes 1 to mu. This is a 4-space
% rotation which is a special case of the general 4-space rotation given
% by Coxeter (see note 2 below).

r = sqrt(unit(mu));
S = r .* S .* r;

end

function Z = beta
% Generate a random value with the beta(a, a) distribution with a = 1.5.
% (1.5 because the parameters for beta are (m - 1)/2 and m = 4.)

% Reference: Ulrich, p161, algorithm SB (see above).
 
while true
    U = 2 .* rand - 1.0; % A uniformly distributed value in (-1, 1).
    V = rand;

    S = U.^2 + V.^2;

    if S <= 1
        Z = 0.5 + U .* V .* sqrt(1 - S) ./ S;
        return
    end
end

end

% Note 1. The call to the internal function beta is equivalent to calling
% the Matlab statistics toolbox function betarnd(1.5,1.5) but we do not use
% betarnd to avoid dependence on the Matlab statistics toolbox, and also
% because what is needed here is a special case that can be very simply
% coded.

% Note 2. The 4-space rotation that maps (1,0,0,0) to mu is derived as
% follows, using Coxeter's Theorem 5.1 and 5.2 (see §5 of Coxeter's paper,
% especially the statement just after Theorem 5.1). The general case of a
% 4-space rotation that rotates q to p is:
%
% p = sqrt(p .* conj(q)) .* q .* sqrt(conj(q) .* p)
%
% This is obtained by taking the square root of the left and right terms of
% Coxeter's formula, in order to make the angle of rotation the same as the
% angle between p and q (whose cosine is the scalar product of p and q).
%
% In the special case used here, q = (1,0,0,0), and the formula reduces to:
%
% p = sqrt(p) .* q .* sqrt(p)
%
% Reference:
%
% H. S. M. Coxeter, "Quaternions and Reflections", American Mathematical
% Monthly, 53(3), pp136-146, March 1946.

% $Id: randvmf.m 1004 2017-11-15 17:14:09Z sangwine $
