function s = field(f, realStr, complexStr, quaternionStr)
% Returns either:
% - 'realStr' if f = 'R'
% - 'complexStr' if field = 'C'  
% - 'quaternionStr' if field = 'H'
%
% Default values for realStr/complexStr/quaternionStr are 'real'/'complex'/'quaternion'
    if nargin < 4
        quaternionStr = 'quaternion';
    end
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
      case 'H'
        s = quaternionStr;
      otherwise
        error(sprintf('Unknown field %s', f));
    end
end
