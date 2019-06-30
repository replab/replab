function str = char(o)
% CHAR Create character array (string).
% (Octonion overloading of standard Matlab function.)

% Note: the Matlab char function converts arrays of numeric values into
% character strings. This is not what this function does, but the Matlab
% guidance on user-defined classes suggests writing a char function and
% a disp/display function. This advice has been followed.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

% There are three cases to be handled. The argument is one of: empty, a
% pure octonion, a full octonion.

if isempty(o)
    str = '[] octonion'; % This case must be handled first because an
    return;                % empty octonion is not scalar, and the
end                        % next check would fail.

if ~isscalar(o)
    error('char cannot handle a vector or a matrix.')
end

f = '%.4g'; % Control over the format of each numeric value.

if ~isempty(s(o.a))
   % There is a scalar part to output:
   t = s(o.a);
   if isreal(t)
       % The scalar part is real. A minus sign will be provided by num2str,
       % if needed.
       
       str = [num2str(t, f) ' '];
   else
       % The scalar part is complex, so just output it in parentheses, with
       % no sign (the real and imaginary parts will have a sign supplied by
       % num2str).
       
       str = ['(' num2str(t, f) ') '];
   end
else
   str = '';
end

% Now output the X/Y/Z/A/B/C/D components, with I/J/K/L/M/N and suitable
% signs and spaces. If the octonion is pure, the X component is not
% preceded by a scalar part, and we treat it differently to avoid a space
% between the sign and the numeric value.

t = x(o.a);
if isempty(s(o.a))
    % No scalar part exists, so the X component is the first piece to be
    % output.
    
    if isreal(t)
        str = [num2str(t, f) ' * I '];
    else
        str = ['(' num2str(t, f) ') * I '];
    end
else
    % There was a preceding scalar part, so we have to space accordingly.
    
    if isreal(t)
        str = [str plusminus(t) ' ' num2str(abs(t), f) ' * I '];
    else
        str = [str '+ (' num2str(t, f) ') * I '];
    end
end

t = y(o.a);
if isreal(t)
    str = [str plusminus(t) ' ' num2str(abs(t), f) ' * J '];
else
    str = [str '+ (' num2str(t, f) ') * J '];
end

t = z(o.a);
if isreal(t)
    str = [str plusminus(t) ' ' num2str(abs(t), f) ' * K '];
else
    str = [str '+ (' num2str(t, f) ') * K '];
end

t = s(o.b);
if isreal(t)
    str = [str plusminus(t) ' ' num2str(abs(t), f) ' * L '];
else
    str = [str '+ (' num2str(t, f) ') * L '];
end

t = x(o.b);
if isreal(t)
    str = [str plusminus(t) ' ' num2str(abs(t), f) ' * M '];
else
    str = [str '+ (' num2str(t, f) ') * M '];
end

t = y(o.b);
if isreal(t)
    str = [str plusminus(t) ' ' num2str(abs(t), f) ' * N '];
else
    str = [str '+ (' num2str(t, f) ') * N '];
end

t = z(o.b);
if isreal(t)
    str = [str plusminus(t) ' ' num2str(abs(t), f) ' * O'];
else
    str = [str '+ (' num2str(t, f) ') * O'];
end
end

function S = plusminus(X)
% Extracts the sign of X and returns the character '-', '+'.

if sign(X) == -1, S = '-'; else S = '+'; end
end

% $Id: char.m 1004 2017-11-15 17:14:09Z sangwine $
