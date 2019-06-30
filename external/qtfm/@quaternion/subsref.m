function b = subsref(a, ss)
% SUBSREF Subscripted reference.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if length(ss) ~= 1
    error('Only one level of subscripting is currently supported for quaternions.');
    % See the notes below under structure indexing.
end

switch ss.type
case '()'
    if length(ss) ~= 1
        error('Multiple levels of subscripting are not supported for quaternions.')
    end
    
    % To implement indexing, we operate separately on the components.
    
    if isempty(a.w)
        xa = a.x; ya = a.y; za = a.z;
        b = a; % Copy the input parameter to avoid calling the constructor.
        b.x = xa(ss.subs{:});
        b.y = ya(ss.subs{:});
        b.z = za(ss.subs{:});
    else
        sa = a.w; xa = a.x; ya = a.y; za = a.z;
        b = a; % Copy the input parameter to avoid calling the constructor.
        b.w = sa(ss.subs{:});
        b.x = xa(ss.subs{:});
        b.y = ya(ss.subs{:});
        b.z = za(ss.subs{:});
    end
case '{}'
    error('Cell array indexing is not valid for quaternions.')
case '.'
    % Structure indexing. TODO.
    %
    % See some notes on this subject in the file subsasgn.m. Here, we would
    % have to support two levels of indexing, such as q.x(1,2) and q(1,2).x
    % which would give the same result. The Matlab help on subsref explains
    % how this would work. We could accept s, x, y, or z as field names and
    % return x(a) etc, but we would really need to pass a second level of
    % subscripting back to subsref recursively. The code below is an
    % interim step to support one level of structure indexing.
    %
    % 15-Sep-2005 Code contributed by T.A. Ell.
    %
    switch ss.subs
        case {'vector', 'v'}    
            b = vee(a);
        case {'scalar', 's', 'w'} % w added as synonym for s, 23 May 2008.
            b = scalar(a); % Calling scalar() rather than s() means an
                           % array of zeros is returned if a is pure.
        case {'x', 'I'}
            b = a.x;
        case {'y', 'J'}
            b = a.y;
        case {'z', 'K'}
            b = a.z;
        case 'imag'
            b = imag(a);
        case {'real'}
            b = real(a);
        otherwise
            error( ['Structure ''.', ss.subs, ''' is not a valid index']);
    end
otherwise
    error('subsref received an invalid subscripting type.')
end

% $Id: subsref.m 1004 2017-11-15 17:14:09Z sangwine $
