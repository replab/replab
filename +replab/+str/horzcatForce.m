function z = horzcatForce(varargin)
% Returns a 'horzcat' of multiple (nargin >= 0) arguments
%
% If any of the arguments is not a horizontal cell array vector,
% reshapes it.
    switch nargin
      case 0
        z = cell(1, 0);
      case 1
        z = varargin{1};
        z = z(:)';
      otherwise
        x = varargin{1};
        x = x(:)';
        rest = varargin(2:end);
        y = replab.str.horzcatForce(varargin{2:end});
        z = horzcat(x, y);
    end
end
