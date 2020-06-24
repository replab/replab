function str = linkHelp(helpCommand, linkText, identifier)
% Returns a HTML link that runs a help command if the console supports HTML, or plain text as a fallback
%
% Args:
%   helpCommand (charstring): Invocation of the help command, possibly including flags, without trailing space
%                             Examples would be 'help -f' or 'help'
%   linkText (charstring): Link text
%   identifier (charstring): Identifier to ask help for
    parts = strsplit(strtrim(helpCommand));
    helpName = parts{1};
    args = {parts{2:end} identifier};
    if replab.globals.consoleUseHTML
        argString = strjoin(cellfun(@replab.infra.quote, args, 'uniform', 0), ',');
        str = sprintf('<a href="matlab: %s(%s)">%s</a>', helpName, argString, linkText);
    else
        str = linkText;
    end
end
