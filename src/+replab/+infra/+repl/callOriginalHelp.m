function [output, docTopic] = callOriginalHelp(coh_argument, coh_printIt, coh_variableName, coh_variableValue)
% Calls the original Matlab/Octave help function.
%
% This function either displays the help content or returns it, as
% instructed by the 'coh_printIt' argument. Printing should not be 
% requested here under matlab/emacs environment.
%
% Args:
%   coh_argument (cell of cell of charstring or object): Argument to call the help function with
%   coh_printIt (boolean): If the help should be displayed, otherwise it is returned
%   coh_variableName (charstring, optional): If given, variable name to insert in the caller workspace
%   coh_variableValue (optional): If given, value corresponding to coh_variableName
%
% Returns
% -------
%   output: charstring
%     The output of the original help command
%   docTopic: charstring
%     Name of the documentation item associated with this topic

    if nargin < 2
        coh_printIt = true;
    end
    if coh_printIt
        assert(~replab.globals.runsInMatlabEmacs);
    end
    
    coh_functionHandle = replab.globals.defaultHelpFunction;
    
    output = '';
    docTopic = '';
    if nargin > 3
        % verify we don't overwrite our own local variables
        assert(~ismember(coh_variableName, {'coh_functionHandle' 'coh_printIt', 'coh_argument' 'coh_variableName' 'coh_variableValue'}));
        eval([coh_variableName ' = coh_variableValue;'])
    end
    if coh_printIt
        coh_functionHandle(coh_argument{1}{:});
    else
        if replab.compat.isOctave
            output = coh_functionHandle(coh_argument{1}{:});
        else
            [output, docTopic] = coh_functionHandle(coh_argument{1}{:});
        end
    end
    
end
