function tell(msg, varargin)
% Debug method to give some verbosity to the decomposition process
    persistent level
    if isempty(level)
        level = 0;
    end
    switch msg
      case 'down'
        level = level + 1;
      case 'up'
        level = level + 1;
      otherwise
        prefix = '';
        for i = 1:level
            prefix = [prefix '..'];
        end
        %disp([prefix sprintf(msg, varargin{:})]);
    end
end
