function m = spherical_mean(q, tol)
% SPHERICAL_MEAN Computes the spherical mean, m, of the elements of q.
% This is a unit quaternion that minimises the sum of geodesic distances
% from m to unit versions of the elements of q. The second (optional)
% parameter is an angular tolerance on the result (in radians).
% q must be a vector (it may be pure or full quaternion(s)) and tol mus be
% a numeric scalar between 0 and 1.

% Copyright (c) 2008, 2016 Stephen J. Sangwine, Nicolas Le Bihan and Salem Said.
% See the file : Copyright.m for further details.

% Reference:
%
% Samuel R. Buss and Jay P. Fillmore, 'Spherical averages and
% applications to spherical splines and interpolation',
% ACM Transactions on Graphics, 20(2), 95-126, April 2001.
% DOI:10.1145/502122.502124.
%
% This paper describes two algorithms for computing the spherical mean (§3),
% and discusses the general principles (in §2). The algorithm used here is
% a modification of Algorithm A1.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 1
    tol = 1e-6; % The tolerance was not specified, supply a default.
                % 1e-6 radians is less than one second of arc.
else
    if (~isnumeric(tol) || ~isscalar(tol))
        error('Second parameter must be a numeric scalar.')
    end
    if tol <= 0 || tol >= 1
        error('Second parameter must be positive and less than 1.')
    end
end

if ~isvector(q)
    % TODO This function could be modified to handle matrices also, cf the
    % Matlab mean function. It would need to work on columns.
    error('First parameter must be a vector.')
end

if ~isreal(q)
    warning('QTFM:information', ...
        'The validity of SPHERICAL_MEAN on complexified quaternions is unknown.')
end

u = unit(q);

% u is now a vector of unit quaternions with the same directions as
% elements of q.

m = mean(u); % Initial value, somewhat arbitrary: take the Euclidean mean.

if abs(abs(m)) < 1e-2 % Use double abs in case m is complexified.
    % The geometric mean is near the centre of the 3/4-sphere.
    warning('Spherical mean may be incorrect (geometric mean is small).')
end

m = unit(m); % Now we can safely normalise the Euclidean mean, having
             % checked that it is non-zero.

P = ispure(q); % We need to know whether q is pure, because if it
               % is, the result must be pure.

while true
    d = (m.^-1) .* u; % d is a vector of quaternions that will postmultiply
                      % m to give u. Each element of d is a geodesic
                      % distance measure since it rotates m to an element
                      % of u. We are aiming to minimise this measure.
    e = exp(mean(log(d))); % This is a correction factor that moves the
                           % current estimate of the spherical mean nearer
                           % to the true value.
    m = m .* e;     % Update the estimate. Note that if q is pure, m should
    if P            % be too, but e will be a full quaternion. Therefore
        m = vee(m); % rounding error can result in a non-zero scalar part
    end             % in m, which we must remove.
    if (abs(angle(e)) <= tol), break, end
end

end

% $Id: spherical_mean.m 1004 2017-11-15 17:14:09Z sangwine $
