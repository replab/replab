function str = linkHelp(helpFunctionName, linkText, helpArg, flags)
% Returns a HTML link that runs a help command if the console supports HTML, or plain text as a fallback
%
% Args:
%   helpFunctionName (charstring): Name of the help command
%   linkText (charstring): Link text
%   helpArg (charstring): Main argument to the help command
%   flags (charstring, row cell vector of charstring, optional): Flag, or flags to pass on to the help command
    if nargin < 3
        lineNumber = [];
    end
    
    if replab.Parameters.consoleUseHTML
        t = linkText;
    else
        t = altText;
    end
    tokens = regexp(t, '(%\w)', 'tokens');
    n = length(tokens);
    args = cell(1, n);
    for i = 1:n
        token = tokens{i};
        token = token{1};
        switch token
          case '%s'
            args{i} = filename;
          case '%d'
            if isempty(lineNumber)
                error('Line number cannot be empty if used in the link message');
            end
            args{i} = lineNumber;
          otherwise
            error('Unknown format specifier');
        end
    end
    t = sprintf(t, args{:});
    if replab.Parameters.consoleUseHTML
        if isempty(lineNumber)
            lineNumber = 1;
        end
        link = sprintf('matlab: opentoline(''%s'', %d)', filename, lineNumber);
        str = ['<a href ="' link '">' t '</a>'];
    else
        str = t;
    end
end
