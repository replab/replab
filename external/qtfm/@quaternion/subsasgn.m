function r = subsasgn(a, ss, b)
% SUBSASGN Subscripted assignment.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2010, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

switch ss.type
case '()'
    if length(ss) ~= 1
        error('Multiple levels of subscripting are not supported for quaternions.')
    end
    
    if ~isa(a, 'quaternion')
        % TODO Consider how Matlab handles r(1) = randn. This works because
        % Matlab creates a real array for the LHS, then assigns to its
        % first element. In the quaternion case, we could handle this by
        % writing a quaternion into r, but would that make sense? Is it
        % better to leave it as is, with an error?
        error('Left side must be a quaternion in subscripted assignment if right side is a quaternion.')    
    end

    if ~isa(b, 'quaternion')

           % Argument a (the left-hand side of the assignment) is a
           % quaternion, but b (the right-hand side) is not (it could be
           % double, for example, zero; or it could be empty). To handle
           % some of these cases, we need to convert b to a quaternion.
           
           if isempty(b)
               % There is a subtle issue here with assigning empty. The
               % following works in double:  A = randn(4,1); A(3) = [].
               % But if you create a variable containing the empty, it
               % won't work. Thus: B = []; A(3) = B; gives the "Subscripted
               % assignment dimension mismatch." error. Therefore to make
               % this work for quaternion arrays, we have to explicitly
               % assign [], rather than attempting to assign the empty that
               % is passed to us in b.
               % TODO There may remain an issue here with the class of the
               % empty array, cf the class of array a.

               if ~isempty(a.w)  
                   sa = a.w; sa(ss.subs{:}) = []; a.w = sa;
               end
               
               xa = a.x; xa(ss.subs{:}) = []; a.x = xa;
               ya = a.y; ya(ss.subs{:}) = []; a.y = ya;
               za = a.z; za(ss.subs{:}) = []; a.z = za;

               r = a;
               return
           end
           
           if isnumeric(b) % but not empty because we dealt with that.
               if isempty(a.w) % a is a pure quaternion ...
                   % ... and we can't handle a numeric in this case, because
                   % we can't make a pure quaternion out of a numeric.
                   error('Cannot convert right-hand (numeric) argument to a pure quaternion.');
               else
                   b = quaternion(b); % Convert b to a quaternion (with implicit
                                      % zero vector part), and carry on.
               end
           else
               error(['Cannot handle right-hand argument of class ', class(b)]);
           end
           
    end
    
    % Both parameters, a and b, are quaternions.
    
    pa = isempty(a.w);
    
    if pa ~= isempty(b.w)
        % TODO Surely we could relax this restriction by converting the
        % left side into a full quaternion by supplying zeros of the
        % appropriate class, and then assigning b. Here is a test case:
        %
        % >>q = randv
        % >>q = q + 0
        %
        % q = 0 - 0.9003 * I + 0.4211 * J - 0.1102 * K   % THIS WORKS OK!
        %
        % >>q = randv(2)
        %
        % q = 2x2 pure quaternion array
        %
        % >>q(1,1) = randq                                 % THIS DOESN'T!
        % Error using quaternion/subsasgn (line 89)
        error(['Left and right sides in subscripted assignment must ', ...
               'both be pure, or both full quaternions.'])
    end
    
    % To implement indexed assignment, we operate separately on the
    % components, replace the components with the modified components,
    % then copy to the output parameter.
     
    if ~pa
        sa = a.w; sa(ss.subs{:}) = b.w; a.w = sa;
    end

    xa = a.x; xa(ss.subs{:}) = b.x; a.x = xa;
    ya = a.y; ya(ss.subs{:}) = b.y; a.y = ya;
    za = a.z; za(ss.subs{:}) = b.z; a.z = za;

    r = a;
case '{}'
    error('Cell array indexing is not valid for quaternions.')
case '.'
    error('Structure indexing is not implemented for quaternions.')
    %
    % Possible use of structure indexing is to implement the following
    % sorts of assignment:
    %
    % q.x = blah
    %
    % However, there are some issues to be considered before implementing
    % such a scheme. Would it work, for example, if q doesn't exist? What
    % should happen if q is pure and q.s = blah occurs? Should q become a
    % full quaternion, or should an error be raised? What about mixed use
    % of array indexing and structure indexing, e.g. q.x(:,1)? Would it
    % work for q.x = q.x .* 2 where the structure indexing occurs on the
    % right hand side as well as on the left. Guess: q.x on the right would
    % be handled by subsref, not subsassgn.
    %
    % (Notes added after discussion with Sebastian Miron, 12 July 2005.)
    %
otherwise
    error('Quaternion subsasgn received an invalid subscripting type.')
end

% $Id: subsasgn.m 1032 2019-04-20 14:20:06Z sangwine $
