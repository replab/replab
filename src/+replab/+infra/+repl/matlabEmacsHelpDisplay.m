function matlabEmacsHelpDisplay(arg, contents)
% Displays a help topic in an Emacs windows when using the matlab-emacs mode
%
% Args:
%   arg (charstring): Argument passed to the help command, will be display in the help frame title
%   contents (charstring): Help contents to display
    fprintf('<EMACSCAP>(*MATLAB Help: %s*)\n%s\n</EMACSCAP>', arg, contents);
end
