function R = overload(F, Q, varargin)
% Private function to implement overloading of Matlab functions. Called to
% apply the function F to the quaternion array Q by operating on components
% of Q with F. F must be a string, giving the name of the function F. The
% calling function can pass this string using mfilename, for simplicity of
% coding. varargin contains optional arguments that are not quaternions.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

H = str2func(F); % A handle to the function designated by F.

R = Q;

if ~isempty(R.w) % Skip the scalar part if it is empty (pure quaternion).
    R.w = H(R.w, varargin{:});
end

R.x = H(R.x, varargin{:});
R.y = H(R.y, varargin{:});
R.z = H(R.z, varargin{:});

% $Id: overload.m 1004 2017-11-15 17:14:09Z sangwine $
