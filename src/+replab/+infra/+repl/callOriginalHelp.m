function output = callOriginalHelp(coh_argument, coh_variableName, coh_variableValue)
% Calls the original Matlab/Octave help function and returns the help contents
%
% Args:
%   coh_argument (charstring or object): Argument to call the help function with
%   coh_variableName (charstring, optional): If given, variable name to insert in the caller workspace
%   coh_variableValue (optional): If given, value corresponding to coh_variableName
%
% Returns:
%   charstring: The output of the original help command
    coh_functionHandle = replab.globals.defaultHelpFunction;
    if nargin > 2
        % verify we don't overwrite our own local variables
        assert(~ismember(coh_variableName, {'coh_functionHandle' 'coh_argument' 'coh_variableName' 'coh_variableValue'}));
        eval([coh_variableName ' = coh_variableValue;'])
        output = coh_functionHandle(coh_argument);
    else
        output = coh_functionHandle(coh_argument);
    end
end
