function o = times(l, r)
% .*  Array multiply.
% (Octonion overloading of standard Matlab function.)
 
% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

ol = isa(l, 'octonion');
or = isa(r, 'octonion');

if ol && or
    
    % Both arguments are octonions.
    
    % Use the Cayley-Dickson formula for the octonion product. This is
    % given in Ward (1997), or the Wikipedia page on octonions as:
    % (a,b)(c,d) = [(ac - d*b), (da + bc*)] where the * denotes a quaternion
    % conjugate.

    o = l; % Create a result, avoiding the use of the constructor.
    o.a = l.a .* r.a - conj(r.b) .* l.b;
    o.b = r.b .* l.a + l.b .* conj(r.a);
else

    % One of the arguments is not an octonion. If it is numeric, we can
    % handle it.
    
    if ol && isa(r, 'numeric')
        o = l; % The left operand is an octonion, so use it for the result
               % to avoid calling the constructor.     
        o.a = o.a .* r; o.b = o.b .* r;
    elseif isa(l, 'numeric') && or
        o = r; % The right operand is an octonion, so use it for the
               % result to avoid calling the constructor.  
        o.a = l .* o.a; o.b = l .* o.b;

    else
        if isa(l, 'quaternion') || isa(r, 'quaternion')
            error('Cannot multiply quaternions and octonions')
        else
            error('Multiplication of an octonion by a non-numeric is not implemented.')
        end
    end
end

% $Id: times.m 1004 2017-11-15 17:14:09Z sangwine $
