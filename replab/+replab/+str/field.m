function s = field(f, realStr, complexStr)
% Returns either:
% - 'realStr' if f = 'R'
% - 'complexStr' if field = 'C'  
%
% Default values for realStr/complexStr are 'real'/'complex'
    if nargin < 3
        complexStr = 'complex';
    end
    if nargin < 2
        realStr = 'real';
    end
    switch f
      case 'R'
        s = realStr;
      case 'C'
        s = complexStr;
      otherwise
        error(sprintf('Unknown field %s', f));
    end
end
