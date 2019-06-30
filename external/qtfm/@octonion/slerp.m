function y = slerp(q1, q2, t)
% SLERP   Spherical linear interpolation.
% Interpolates between two octonions q1 and q2, using parameter t to
% determine how far along the 'arc' to position the result (interior to the
% arc or exterior). t = 0 corresponds to q1, t = 1 corresponds to q2, and t
% = 0.5 will give a value mid-way between q1 and q2. q1, q2 and t must have
% the same size, unless one or more is scalar.

% Copyright (c) 2008, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% The third parameter, t, gives the 'distance' along the 'arc' between the
% octonions, 0 representing q1 and 1 representing q2. If q1 and q2 are
% unit pure octonions, the interpolation is along a great circle of the
% sphere between the points represented by q1 and q2. If q1 and q2 are unit
% full octonions, the interpolation is along the 'arc' on the 8-sphere:
% this means the result is an octonion which represents a rotation
% intermediate between the two rotations represented by q1 and q2. If the
% first two parameters are not unit octonions, then there is also
% interpolation in modulus.

% q1 and q2 are not restricted to be real octonions. Since a unit complex
% octonion has a real modulus, the check that the first two parameters
% have unit modulus works as coded.

% Interpretation.
%
% The slerp function can be simply understood in terms of the ratio of two
% vectors (pure octonions). The ratio is the octonion that rotates one
% vector into the other. Taking a fractional power of this rotation and
% then multiplying it by the first vector obviously gives a vector which is
% part way along the arc between the two octonions. The ratio is computed
% using the multiplicative inverse. If the two octonions are full, then
% their ratio again gives an octonion which multiplies one to give the
% other, but this time in 8-space, including, if the moduli are not unity,
% the scale factor needed to scale one into the other.

% Reference:
%
% Ken Shoemake, 'Animating rotation with quaternion curves', SIGGRAPH
% Computer Graphics, 19 (3), July 1985, 245-254, ACM, New York, USA.
% DOI:10.1145/325165.325242.

narginchk(3, 3), nargoutchk(0, 1)

if ~isnumeric(t) || ~isreal(t)
    error('Third parameter must be real and numeric.');
end

if ~(all(size(q1) == size(q2)) || isscalar(q1) || isscalar(q2))
    error(['First two parameters cannot be of different sizes unless' ...
          ' one is a scalar.']);
end

if ~isscalar(t)
    
    % t is an array. This is only possible if it matches the size of one of
    % the first two parameters, or both of the first two parameters are
    % scalars.
    
    % TODO If t has the same size as one of the other two, then each value
    % in t will be used to raise one value to a power. This may not be what
    % the user expects, if say, q1 is scalar, and q2 is a vector the same
    % size as t. Maybe needs a little more care. Also, t could have one
    % more dimension than q1, q2, with the extra dimension being the
    % interpolating 'vector'. Again needs a bit more thought.
    
    if ~(all(size(q1) == size(t)) || all(size(q2) == size(t)) || ...
         (isscalar(q1) && isscalar(q2)) ...
        )
        error(['Third parameter cannot be an array unless' ...
               ' the first two are scalars, or it has the'...
               ' same size as one of the first two parameters.']);
    end
end

r12 = q1.^-1 .* q2; % This is the octonion that turns q1 into q2.

if any(any(t < 0)) || any(any(t > 1))

    % Up to February 2016 there was a restriction here to the range 0 to 1,
    % raising an error if any element of t fell outside this range. However
    % the restriction is not necessary: a value greater than 1 interpolates
    % outside the arc between q1 and q2, and similarly for a value less
    % than zero. But it does not make sense to use arbitrarily large
    % values, so we check that the values fall within half a revolution to
    % within an arbitrary rounding error.

    if any(any(abs(angle(r12) .* t) - pi > 5 .* eps))
        error('Slerp angle is outside the sensible range of (-pi, pi).');
    end
end

y = q1 .* (r12 .^ t);

if ispure(q1) && ispure(q2) % If both input parameters were pure, so must
    y = vee(y);             % be the result, by definition, so remove any
end                         % scalar part.

end

% $Id: slerp.m 1004 2017-11-15 17:14:09Z sangwine $
