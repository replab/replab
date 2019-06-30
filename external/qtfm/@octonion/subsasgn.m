function r = subsasgn(a, ss, b)
% SUBSASGN Subscripted assignment.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

switch ss.type
case '()'
    if length(ss) ~= 1
        error('Multiple levels of subscripting are not supported for octonions.')
    end
    
    if ~isa(a, 'octonion')
        error('Left side must be an octonion in subscripted assignment if right side is an octonion.')    
    end

    if ~isa(b, 'octonion')

        if isempty(b)
            % There is a subtle issue here with assigning empty. The
            % following works in double:  A = randn(4,1); A(3) = [].
            % But if you create a variable containing the empty, it
            % won't work. Thus: B = []; A(3) = B; gives the "Subscripted
            % assignment dimension mismatch." error. Therefore to make
            % this work for octonion arrays, we have to explicitly
            % assign [], rather than attempting to assign the empty that
            % is passed to us in b.
            % TODO There may remain an issue here with the class of the
            % empty array, cf the class of array a.
            
            r = a; % Copy the input array to create an octonion result.
            
            r.a = subsasgn(a.a, ss, []); % These calls are to the quaternion
            r.b = subsasgn(a.b, ss, []); % subsagn function, of course.
            return
        end
           
        % Argument a (the left-hand side of the assignment) is an octonion,
        % but b (the right-hand side) is not (it could be double, for example,
        % zero). To handle this case, we need to convert b to an octonion.
        
        if isnumeric(b) % but not empty because we dealt with that.
            if ispure(a) % a is a pure octonion ...
                % ... and we can't handle a numeric in this case, because
                % we can't make a pure octonion out of a numeric.
                error('Cannot convert right-hand (numeric) argument to a pure octonion.');
            else
                b = octonion(b); % Convert b to an octonion (with implicit
                                 % zero vector part), and carry on.
            end
        else
            error(['Cannot handle right-hand argument of class ', class(b)]);
        end

    end
    
    pa = ispure(a);
    if pa ~= ispure(b)
        error('Left and right sides in subscripted assignment must both be pure, or both full octonions.')
    end
    
    % To perform subscripted assignment, we split the octonion into two
    % quaternions and operate on the two quaternions separately. Finally we
    % re-assemble the two quaternion results to make the octonion result.
    
    r = a; % Copy the input array to create an octonion result.
    r.a = subsasgn(a.a, ss, b.a); % These calls are to the quaternion
    r.b = subsasgn(a.b, ss, b.b); % subsagn function, of course.
case '{}'
    error('Cell array indexing is not valid for octonions.')
case '.'
    error('Structure indexing is not implemented for octonions.')
otherwise
    error('Octonion subsasgn received an invalid subscripting type.')
end

% $Id: subsasgn.m 1004 2017-11-15 17:14:09Z sangwine $
