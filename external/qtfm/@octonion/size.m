function varargout = size(o, dim)
% SIZE   Size of matrix.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

switch nargout
    case 0
        switch nargin
            case 1
                size(o.a)
            case 2
                size(o.a, dim)
            otherwise
                error('Incorrect number of input arguments.')
        end
    case 1
        switch nargin
            case 1
                varargout{1} = size(o.a);
            case 2
                varargout{1} = size(o.a, dim);
            otherwise
                error('Incorrect number of input arguments.')
        end
    case 2
        switch nargin
            case 1
                [varargout{1}, varargout{2}] = size(o.a);
            case 2
                error('Unknown command option.'); % Note 1.
            otherwise
                error('Incorrect number of input arguments.')
        end
    otherwise
        switch nargin
            case 1
                d = size(o.a);         
                for k = 1:length(d)
                    varargout{k} = d(k); %#ok<AGROW>
                end
            case 2
                error('Unknown command option.'); % Note 1.                
            otherwise
                error('Incorrect number of input arguments.')
        end
end

% Note 1. Size does not support the calling profile [r, c] = size(q, dim),
% or [d1, d2, d3, .... dn] = size(q, dim). The error raised is the same as
% that raised by the built-in Matlab function.

% $Id: size.m 1004 2017-11-15 17:14:09Z sangwine $
