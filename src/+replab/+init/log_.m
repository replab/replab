function log_(level, message, varargin)
    if level <= replab.globals.verboseInit
        fprintf([message '\n'], varargin{:});
    end
end
