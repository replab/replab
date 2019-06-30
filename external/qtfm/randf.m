function S = randf(mu, kappa, varargin)
% RANDF   Returns unit vectors (pure quaternions) distributed on the sphere
% according to a Fisher distribution (a special case of the more general
% von Mises-Fisher distributions). See also RANDVMF for the 4-sphere case.
%
% The first parameter must be a pure quaternion, defining the direction of
% the centre of the distribution. The second parameter is the concentration
% parameter which controls the spread of the distribution on the sphere. It
% must be non-negative. A value of zero results in a uniform distribution
% on the sphere. Larger values result in greater concentration of the
% distribution in the mean direction mu.
%
% The remaining parameters are as for the Matlab function rand (q.v.).  The
% result may be scalar, vector, matrix or array depending on the parameters
% supplied.  Each pure quaternion returned is the result of at least two
% calls on rand, and two calls on randn, and hence randf modifies the state
% of the generator used by both rand and randn. To initialise the generator
% or control the choice of generator, use rand and/or randn.

% This code was rearranged in January 2016 (CVS revision 1.7) to vectorise
% the calculation of S. A side effect of this change is that calls to rand
% and randn are done in a different order to previous versions, which means
% that the function will not produce the same random values as it did, even
% if the random number generator has the same initial value.

% Copyright (c) 2008, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% References:
%
% [The paper by Fisher is the definitive paper on the distribution in
% 3-dimensions. Mardia and Jupp discuss it, but refer to Wood for details
% of how to generate a set of samples from the distribution. Wood's work is
% based on the paper by Ulrich. The report by Dhillon and Sra was used to
% code the implementation below, but they gave an algorithm for the general
% m-dimensional case. We have simplified the code for the case of three
% dimensions and also made a correction -- their Figure 4 does not make use
% of c after computing it. Wood's paper reveals that it should be used in
% the test on Z and U.]
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

if ~isa(mu, 'quaternion') || ~ispure(mu)
    error('The first parameter (mean direction) must be a pure quaternion.')
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

% TODO Each element of w requires a small number of iterations, possibly
% only one. This makes it hard to vectorise the code, and for the moment,
% we adopt a non-vectorised approach. However, the part of the code that
% does the main computation is now vectorised and to vectorise the
% remaining code would not be too difficult. Profiling reveals that very
% little time is spent in this loop (check).

for i = 1:numel(w)
    while true
        z = rand; % This is in fact betarnd(1,1), which is the same as rand
        u = rand; % because we are working in three dimensions.
        t = (1 - (1 + b) .* z) ./ (1 - (1 - b) .* z);
        if kappa * t + 2 * log(1 - x0 .* t) - c >= log(u), break, end
    end
    w(i) = t;
end

% Now compute S in one vectorised step, using the values of w.

% We now need to create a point on the unit circle, uniformly
% distributed. This can be done with two Gaussian random values,
% normalised as if they were a complex number.

v = sqrt(1 - w.^2) .* sign(complex(randn(size(w)), randn(size(w))));
S = quaternion(real(v), imag(v), w);

% S now has elements with the Fisher distribution, but mean direction
% (0,0,1) i.e. aligned with qk or the z-axis. We must now rotate the
% elements of S by the rotation that takes qk to mu. NB mu is not
% necessarily a unit pure quaternion, therefore the rotation is coded using
% the inverse on the right and not the conjugate.

q = sqrt(mu ./ qk);
S = vector(q .* S .* inv(q)); % We take the vector part, because the
                              % rotation gives a zero scalar part
                              % (theoretically).

% $Id: randf.m 1004 2017-11-15 17:14:09Z sangwine $
