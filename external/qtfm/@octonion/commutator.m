function c = commutator(a, b, option)
% COMMUTATOR Computes the commutator of the first two parameters. The
% option parameter selects between two definitions (and meanings) of the
% term:
%
% 'diff' means the pure octonion, c = ab - ba;
% 'prod' means the unit octonion, c, that satisfies: abc = ba.
%
% The default, if no option is specified, is 'diff'.

% Copyright (c) 2013, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% References:
%
% David Nelson, editor, The Penguin Dictionary of Mathematics. Penguin 
% Books, London, third edition, 2003. 'Commutator'
%
% Commutator, Encyclopedia of Mathematics. Accessed 20 August 2015.
% http://www.encyclopediaofmath.org/index.php?title=Commutator&oldid=31490

% Note: there are two definitions of the commutator. One is the difference
% between the product ab and ba, the other is an octonion that multiplies
% ab to give ba. Both are implemented here.

narginchk(2, 3), nargoutchk(0, 1)

if nargin == 2, option = 'diff'; end % Supply the default option.

if ~isa(a, 'octonion') || ~isa(b, 'octonion')
    error('First two parameters must be octonions')
end

if ~ischar(option)
    error('Third parameter must be a character string')
end

if strcmp(option, 'diff')
    
    % This is the most common definition of the commutator (see references
    % above). By definition this is a pure octonion, so we take the vector
    % part to eliminate the scalar part which will be zero to within
    % rounding error.
    
    c = vector(a .* b - b .* a);
    
elseif strcmp(option, 'prod')
    
    % This is the unit octonion that multiplies the product of a and b to
    % give the product of b and a. The result could be calculated as:
    % a.^-1 .* b.^-1 .* a .* b, but taking the inverse requires division by
    % the norm, and it is faster to normalise the result below obtained
    % using conjugates, thus performing the norm and division only once.
    
    c = unit(conj(b) .* conj(a) .* b .* a);  
else
    error('Third parameter must be either ''diff'' or ''prod''')
end

end

% $Id: commutator.m 1004 2017-11-15 17:14:09Z sangwine $
