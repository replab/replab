function res = trileanOr(varargin)
% Returns the 'or' of trilean values
%
% See https://en.wikipedia.org/wiki/Three-valued_logic
%
% True and false are represented by `logical` values, while unknown
% is represented by ``[]``.
%
% Args:
%   varargin: Trilean values
%
% Returns:
%   true or false or []: Result of the 'or' operation
    if nargin == 0
        res = false; % neutral element for or
        return
    end
    isUnknown = cellfun(@(x) isempty(x), varargin);
    if any(isUnknown)
        res = [];
        return
    end
    res = any(cell2mat(varargin));
end
