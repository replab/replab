function s = realComplex(field, realStr, complexStr)
% Returns either 'realStr' if field = 'R' or 'complexStr' if field = 'C'  
%
% Default values for realStr/complexStr are 'real'/'complex'
    if nargin < 3
        complexStr = 'complex';
    end
    if nargin < 2
        realStr = 'real';
    end
    switch field
      case 'R'
        s = realStr;
      case 'C'
        s = complexStr;
      otherwise
        error(sprintf('Unknown field %s', field));
    end
end
