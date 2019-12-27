function str = linkOpen(linkText, altText, filename, lineNumber)
% Returns a HTML link that opens a source file in the editor if the console supports HTML, or alternative text
%
% Args:
%   linkText (charstring): Link text
%                         In there, ``%s`` is replaced by `filename` and ``%d`` by `lineNumber`
%   altText (charstring): Alternative text if the console output does not support HTML
%                         In there, ``%s`` is replaced by `filename` and ``%d`` by `lineNumber`
%   filename (charstring): Full path of the filename to open
%   lineNumber (integer, optional): Line number to open the file at
    if nargin < 3
        lineNumber = [];
    end
    
    if replab.settings.consoleUseHTML
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
    if replab.settings.consoleUseHTML
        if isempty(lineNumber)
            lineNumber = 1;
        end
        link = sprintf('matlab: opentoline(''%s'', %d)', filename, lineNumber);
        str = ['<a href ="' link '">' t '</a>'];
    else
        str = t;
    end
end
