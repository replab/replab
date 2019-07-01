function z = uniqueNames(varargin)
% Merges cell arrays of strings given as arguments (nargin >= 0)
%
% Makes sure that the returned cell array is a row cell array,
% and removes duplicates
    z = replab.str.horzcatForce(varargin{:});
    z = unique(z);
end
