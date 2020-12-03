function log(level, message, varargin)
% Displays a diagnostic message on the console output
%
% Additional arguments can be function handles, to be evaluated only when the message is actually displayed.
%
% Args:
%   level (integer): Verbosity level beginning at which to show the message (see `+replab.+globals.verbosity`)
%   message (charstring): Message to display, can include sprintf-style formatting directives
    if level <= replab.globals.verbosity
        args = varargin;
        for i = 1:length(args)
            f = args{i};
            if isa(f, 'function_handle')
                args{i} = f();
            end
        end
        message = sprintf(message, args{:});
        disp(message);
    end
end
