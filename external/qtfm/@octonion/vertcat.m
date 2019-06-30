function c = vertcat(varargin)
% VERTCAT Vertical concatenation.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

c = varargin{1};

if length(varargin) == 1
    % We have been passed one argument, nothing to do, except check that it
    assert(isa(c, 'octonion'), ... % really is an octonion.
        'Octonion vertcat called with single non-octonion argument!');
    return
else
    if ~isa(c, 'octonion')
        c = octonion(c); % If the first argument is not already an octonion
                         % we need to make it into one, if we can (if we
                         % can't an error will be raised by the constructor,
                         % which will do for now).
    end
end

for i = 2:length(varargin)
    d = varargin{i};
    if isempty(d), continue, end % We can skip an empty argument.
    if isa(d, 'numeric')
       d = octonion(d); % Promote d to an octonion. Doing this first
                        % simplifies the next step.
    end
    if isa(d, 'octonion')
        c.a = [c.a; d.a]; % These calls are to the quaternion
        c.b = [c.b; d.b]; % vertical concatenate.
    else
       cd = class(d);
       error(['Cannot concatenate ', cd, ' with octonion.']); 
    end
end
end

% $Id: vertcat.m 1004 2017-11-15 17:14:09Z sangwine $
