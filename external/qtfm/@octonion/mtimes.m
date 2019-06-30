function o = mtimes(l, r)
% *  Matrix multiply.
% (Octonion overloading of standard Matlab function.)
 
% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

ol = isa(l, 'octonion');
or = isa(r, 'octonion');

% TODO This code avoids using the constructor (for good reason - to avoid
% time-costly error checks which we don't need) but it makes the mistake of
% doing this by copying one of the input parameters. Since this is a matrix
% multiplication, this can result in an initial value for o which is of
% incorrect size. Once assignments are made to o this will be fixed, but
% the change of size may also be costly. This requires study and a careful
% decision made whether to change it. Note that calling a zeros function
% would also call the constructor, which would also be costly in time. Any
% proposed change will require execution time checks on the test code since
% this is such a core function within the toolbox. NB The same
% considerations do not apply to the elementwise product in times.m, since
% in that case the input parameters and the result must have the same size.

if ol && or
    
    % Both arguments are octonions.
    
    % Use the Cayley-Dickson formula for the octonion product. This is
    % given in Ward (1997), or the Wikipedia page on octonions as:
    % (a,b)(c,d) = [(ac - d*b), (da + bc*)] where the * denotes a quaternion
    % conjugate.

    o = l; % Create a result, avoiding the use of the constructor.

    % The formula quoted above is correct for single octonions, but it
    % won't work for matrices because the dimensions may be incompatible.
    % l and r must be conformable, not necessarily the same size, and
    % therefore we must multiply l on the left every time. The solution is
    % to use conjugation, as shown below, to re-order the quaternion
    % products.
    
    o.a = l.a * r.a                   - conj(conj(l.b) * r.b);
    o.b = conj(conj(l.a) * conj(r.b)) +       l.b * conj(r.a);

else

    % One of the arguments is not an octonion. If it is numeric, we can
    % handle it.
    
    if ol && isa(r, 'numeric')
        o = l; % The left operand is an octonion, so use it for the result
               % to avoid calling the constructor.     
        o.a = o.a * r; o.b = o.b * r;
    elseif isa(l, 'numeric') && or
        o = r; % The right operand is an octonion, so use it for the
               % result to avoid calling the constructor.  
        o.a = l * o.a; o.b = l * o.b;

    else
        if isa(l, 'quaternion') || isa(r, 'quaternion')
            error('Cannot multiply quaternions and octonions')
        else
            error('Multiplication of an octonion by a non-numeric is not implemented.')
        end
    end
end

% $Id: mtimes.m 1004 2017-11-15 17:14:09Z sangwine $
