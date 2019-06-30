function S = randvo(varargin)
% RANDVO   Creates uniformly distributed unit pure octonions (points on the
% 7-sphere). Accepts the same parameters as the Matlab function rand
% (q.v.). The result may be scalar, vector, matrix or array depending on
% the parameters supplied. Each octonion returned is the result of seven
% calls on rand, and hence randvo modifies the state of the generator used
% by rand. To initialise the generator or control the choice of generator,
% use rand.

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO Consider a better algorithm, and even whether this algorithm does
% what it says on the tin. The corresponding function for quaternions uses
% a faster algorithm due to James Arvo, which may be adaptable for
% octonions.

S = octonion(randn(varargin{:}), randn(varargin{:}), ...
             randn(varargin{:}), randn(varargin{:}), ...
             randn(varargin{:}), randn(varargin{:}), ...
             randn(varargin{:}));
S = unit(S);

end

% $Id: randvo.m 1004 2017-11-15 17:14:09Z sangwine $
