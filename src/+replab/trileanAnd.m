function res = trileanAnd(varargin)
% Returns the 'and' of trilean values
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
%   true or false or []: Result of the 'and' operation
    if nargin == 0
        res = true; % neutral element for and
        return
    end
    isUnknown = cellfun(@(x) isempty(x), varargin);
    if any(isUnknown)
        res = [];
        return
    end
    res = all(cell2mat(varargin));
end
