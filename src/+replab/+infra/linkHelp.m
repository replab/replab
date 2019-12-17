function str = linkHelp(helpFunctionName, linkText, helpArg, flags)
% Returns a HTML link that runs a help command if the console supports HTML, or plain text as a fallback
%
% Args:
%   helpFunctionName (charstring): Name of the help command
%   linkText (charstring): Link text
%   helpArg (charstring): Main argument to the help command
%   flags (charstring, row cell vector of charstring, optional): Flag, or flags to pass on to the help command
    if nargin < 4
        flags = {};
    end
    if isa(flags, 'char')
        flags = {flags};
    end
    if replab.Parametres.consoleUseHTML
        args = horzcat(flags, {helpArg});
        args = strjoin(cellfun(@(a) ['''' a ''''], args, 'uniform', 0), ',');
        t = sprintf('<a href="matlab: %s(%s)>%s</a>', helpFunctionName, args, linkText);
    else
        t = linkText;
    end
end
