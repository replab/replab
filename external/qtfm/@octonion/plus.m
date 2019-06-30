function o = plus(l, r)
% +   Plus.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Three cases have to be handled:
%
% l is an octonion, r is not,
% r is an octonion, l is not,
% l and r are both octonions.

% An additional complication is that the parameters may be empty (numeric,
% one only) or empty (octonion, one or both). Matlab adds empty + scalar
% to give empty: [] + 2 gives [], not 2, but raises an error (Matrix
% dimensions must agree) on attempting to add empty to an array. This
% behaviour is copied here, whether the empty parameter is an octonion or
% a numeric empty. The result is always an octonion empty.

ol = isa(l, 'octonion');
or = isa(r, 'octonion');

ls = size(l);
rs = size(r);

% Now check for the case where one or other argument is empty and the other
% is scalar. In this case we choose to return an empty octonion, similar
% to the behaviour of Matlab's plus function.

sl = isscalar(l);
sr = isscalar(r);
el = isempty(l);
er = isempty(r);

if (el && sr) || (sl && er)
    o = octonion.empty;
    return
end

% In the quaternion function of the same name, the class of the components
% of the left and right parameters is checked. We omit this here, for
% simplicity, leaving it to the quaternion function to detect this and
% raise an error.

% Having now eliminated the cases where one parameter is empty and the
% other is scalar, the parameters must now either match in size, or one
% must be scalar. The concept of 'matching size' changed with Matlab R2016b
% as from that release Matlab accepts singleton dimensions as matching
% non-singleton dimensions (implicit singleton expansion). Previously there
% were checks here on the sizes, but the complexity of doing this with
% implicit expansion made it simpler to remove the checks and allow Matlab
% to perform them when the built-in plus function is called.
%
% It could be that both parameters are empty arrays of size [0, n] or
% [n, 0] etc, in which case we must return an empty array of the same size
% to match Matlab's behaviour.

if el && er
   % We must return an empty (octonion) matrix.
   if ol && or
       o = l; % r would do just as well, since both are octonions.
   elseif ol
       o = l; % This must be l because r isn't an octonion.
   else
       o = r; % This must be r because l isn't an octonion.
   end
   return
end

if ol && or

    o = l;
    
    o.a = l.a + r.a;
    o.b = l.b + r.b;
    
elseif isa(r, 'numeric')
    
    % The left parameter must be an octonion, otherwise this function
    % would not have been called.
    
    o = l + octonion(r);
    
elseif isa(l, 'numeric')
    
    % The right parameter must be an octonion, otherwise this function
    % would not have been called.
    
    o = octonion(l) + r;
    
else
    if isa(l, 'quaternion') || isa(r, 'quaternion')
        error('Cannot add quaternions and octonions')
    else
        error('Unhandled parameter types in octonion function +/plus')
    end
end

% $Id: plus.m 1004 2017-11-15 17:14:09Z sangwine $
