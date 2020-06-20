function output = callOriginalHelp(coh_function_handle, coh_argument, coh_variableName, coh_variableValue)
% Calls the original Matlab/Octave help function and returns the help contents
%
% Args:
%   coh_function_handle (function_handle): Handle to the original help function taking a single argument
%   coh_argument (charstring or object): Argument to call the help function with
%   coh_variableName (charstring, optional): If given, variable name to insert in the caller workspace
%   coh_variableValue (optional): If given, value corresponding to coh_variableName
%
% Returns:
%   charstring: The output of the original help command
    if nargin > 2
        assert(~ismember(coh_variableName, {'coh_function_handle' ' coh_argument' 'coh_variableName' 'coh_variableValue'}));
        eval([coh_variableName ' = coh_variableValue;'])
        output = coh_function_handle(coh_argument);
    else
        output = coh_function_handle(coh_argument);
    end
end
