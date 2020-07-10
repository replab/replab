function [output, docTopic] = callOriginalHelp(coh_argument, coh_variableName, coh_variableValue)
% Calls the original Matlab/Octave help function and returns the help contents
%
% Args:
%   coh_argument (cell of cell of charstring or object): Argument to call the help function with
%   coh_variableName (charstring, optional): If given, variable name to insert in the caller workspace
%   coh_variableValue (optional): If given, value corresponding to coh_variableName
%
% Returns
% -------
%   output: charstring
%     The output of the original help command
%   docTopic: charstring
%     Name of the documentation item associated with this topic

    coh_functionHandle = replab.globals.defaultHelpFunction;
    if nargin > 2
        % verify we don't overwrite our own local variables
        assert(~ismember(coh_variableName, {'coh_functionHandle' 'coh_argument' 'coh_variableName' 'coh_variableValue'}));
        eval([coh_variableName ' = coh_variableValue;'])
        [output, docTopic] = coh_functionHandle(coh_argument{1}{:});
    else
        [output, docTopic] = coh_functionHandle(coh_argument{1}{:});
    end
end
