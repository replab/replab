function d = associator(a, b, c, option)
% ASSOCIATOR Computes the associator of the three octonions a, b, c. The
% option parameter selects between two definitions (and meanings) of the
% term:
%
% 'diff' means the pure octonion, d = (ab)c - a(bc);
% 'prod' means the octonion, d, that satisfies: ((ab)c)d = a(bc).
%
% The default, if no option is specified, is 'diff'.

% References:
%
% Richard D. Schafer, 'An Introduction to Non-Associative Algebras',
% Academic Press, 1966. Page 13.
%
% Sangwine, S. J., ?Octonion associators?, e-print
% http://arxiv.org/abs/1509.07718, 25 September 2015. 

% Copyright (c) 2013, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(3, 4), nargoutchk(0, 1)

if nargin == 3, option = 'diff'; end % Supply the default option.

if ~isa(a, 'octonion') || ~isa(b, 'octonion') || ~isa(c, 'octonion')
    error('First three parameters must be octonions')
end

if ~ischar(option)
    error('Fourth parameter must be a character string')
end

if strcmp(option, 'diff')
    
    % This is the definition of the associator given by Schafer (see
    % references above). By definition this is a pure octonion, so we take
    % the vector part to eliminate the scalar part which will be zero to
    % within rounding error.
    
    d = vector((a .* b) .* c - a .* (b .* c));
    
elseif strcmp(option, 'prod')
    
    % This is the definition of the associator applicable to quasi-groups
    % and therefore octonions See Sangwine - reference above. It is the
    % unit octonion that multiplies the product of a, b and c evaluated in
    % one order, to give the product evaluated in the other order.
    
    d = ((a .* b) .* c).^-1 .* (a .* (b .* c));
else
    error('Fourth parameter must be either ''diff'' or ''prod''')
end

end

% $Id: associator.m 1004 2017-11-15 17:14:09Z sangwine $
