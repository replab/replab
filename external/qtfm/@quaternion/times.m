function q = times(l, r)
% .*  Array multiply.
% (Quaternion overloading of standard Matlab function.)
 
% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

ql = isa(l, 'quaternion');
qr = isa(r, 'quaternion');

if ql && qr
    
    % Both arguments are quaternions. There are four cases to handle,
    % dependent on whether the arguments are pure or full, the first three
    % being subsets of the full quaternion product.

    pl = isempty(l.w);
    pr = isempty(r.w);

    if pl
        if pr
            % Both arguments are pure quaternions.
            q = l; % Create a result by copying one of the arguments. This
                   % avoids a costly call to the constructor.
            q.w = -(l.x .* r.x + l.y .* r.y + l.z .* r.z);
            q.x =                l.y .* r.z - l.z .* r.y;
            q.y = - l.x .* r.z              + l.z .* r.x;
            q.z =   l.x .* r.y - l.y .* r.x;
        else
            % The right argument is full, but the left is pure.
            q = r;
            q.w = -(l.x .* r.x + l.y .* r.y + l.z .* r.z);
            q.x =   l.x .* r.w + l.y .* r.z - l.z .* r.y;
            q.y = - l.x .* r.z + l.y .* r.w + l.z .* r.x;
            q.z =   l.x .* r.y - l.y .* r.x + l.z .* r.w;
        end
    else
        if pr
            % The left argument is full, but the right is pure.
            q = l;
            q.w = - l.x .* r.x - l.y .* r.y - l.z .* r.z;
            q.x =   l.w .* r.x              + l.y .* r.z - l.z .* r.y;
            q.y =   l.w .* r.y - l.x .* r.z              + l.z .* r.x;
            q.z =   l.w .* r.z + l.x .* r.y - l.y .* r.x;
        else
            % Both arguments are full quaternions. The full monty.
            q = l;
            % TODO Profiling of the test code reveals that the line below
            % takes twice as long to execute as the other three lines. Not
            % clear why, removing parentheses saved only a fraction of the
            % extra time.
            q.w =  l.w .* r.w - l.x .* r.x - l.y .* r.y - l.z .* r.z;
            q.x =  l.w .* r.x + l.x .* r.w + l.y .* r.z - l.z .* r.y;
            q.y =  l.w .* r.y - l.x .* r.z + l.y .* r.w + l.z .* r.x;
            q.z =  l.w .* r.z + l.x .* r.y - l.y .* r.x + l.z .* r.w;
        end
    end
   
else

    % One of the arguments is not a quaternion. If it is numeric or logical
    % we can handle it.
    
    if ql && (isa(r, 'numeric') || isa(r, 'logical'))
        q = l; % The left operand is a quaternion, so use it for the result
               % to avoid calling the constructor.     
        if ~isempty(l.w), q.w = q.w .* r; end
        q.x = q.x .* r; q.y = q.y .* r; q.z = q.z .* r;
    elseif (isa(l, 'numeric') || isa(l, 'logical')) && qr
        q = r; % The right operand is a quaternion, so use it for the
               % result to avoid calling the constructor.  
        if ~isempty(r.w), q.w = l .* q.w; end
        q.x = l .* q.x; q.y = l .* q.y; q.z = l .* q.z;

    else
        error('Multiplication of a quaternion by a non-numeric is not implemented.')
    end
end

% $Id: times.m 1033 2019-04-20 14:35:56Z sangwine $
