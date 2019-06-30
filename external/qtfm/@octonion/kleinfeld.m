function k = kleinfeld(w, x, y, z)
% KLEINFELD  Computes the Kleinfeld product of four octonions.

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(4, 4), nargoutchk(0, 1)

% All four parameters must be octonions. This could possibly be relaxed to
% permit reals, but for the moment, until the utility of this function is
% apparent, we adopt the simplest possible error checking.

if ~isa(w, 'octonion') || ~isa(x, 'octonion') || ...
   ~isa(y, 'octonion') || ~isa(z, 'octonion')
    error('All four parameters must be octonion')
end

% We don't check sizes, but maybe we should....

% Reference:
%
% Erwin Kleinfeld, 'A characterization of the Cayley numbers', pp126-143 in
% A. A. Albert (ed), 'Studies in Modern Algebra', Volume 2, Prentice-Hall
% for The Mathematical Association of America, 1963.

k =        associator(w .* x, y, z) ...
    - x .* associator(w, y, z) ...
    -      associator(x, y, z) .* w;

end

% $Id: kleinfeld.m 1004 2017-11-15 17:14:09Z sangwine $