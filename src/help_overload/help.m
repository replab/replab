function help(varargin)
% Displays help on a function/object/...
%
% Args:
%   varargin: - if a single element, provides help for this element,
%               whether it is a variable, object, function, class,
%               package, ...
%             - if two elements, with one of them being the string '-f'
%               or '-full', provides full help on this element
%
% Note:
%   This function overloads matlab's 'help' function. Renaming this file
%   does not alter its functionality.
%
% Note:
%   The display of RepLAB objects is driven by the liquid-like templates
%   in ``src/+replab/+infra``. Those templates support Sphinx-like reference
%   syntax using backticks.

    contents = {};

    if nargin == 1 && isequal(varargin{1}, '--clear')
        replab.globals.codeBase([]);
        return
    end

    % Are we in full mode or not?
    fullMode = false;
    if nargin == 1
        arg = varargin{1};
    elseif nargin == 2
        if isequal(varargin{1}, '-f') || isequal(varargin{1}, '--full')
            fullMode = true;
            arg = varargin{2};
        elseif isequal(varargin{2}, '-f') || isequal(varargin{2}, '--full')
            fullMode = true;
            arg = varargin{1};
        end
    else
        error('Too many arguments');
    end

    % extract own name
    [~, helpFunctionName, ~] = fileparts(which(mfilename));
    if ischar(arg)
        parts = strsplit(arg, '.');
        isVar = 0;
        try
            varName = parts{1};
            isVar = evalin('caller', ['exist(''' varName ''')']);
            if isVar == 1
                varType = evalin('caller', ['class(' varName ')']);
                varValue = evalin('caller', varName);
            end
        catch
            isVar = 0;
            varName = [];
        end
        if isVar == 1
            % we were passed a string, and the first part is a variable name in the caller workspace
            if replab.compat.startsWith(varType, 'replab')
                % it is a RepLAB object
                contents = replab.infra.repl.lookupHelp(helpFunctionName, fullMode, strjoin({varType parts{2:end}}, '.'), varName);
            else
                contents = replab.infra.repl.callOriginalHelp(arg, varName, varValue);
            end
        else
            if replab.compat.startsWith(arg, 'replab') && ~isequal(arg, 'replab_init') && ~isequal(arg, 'replab_easyinstall')
                contents = replab.infra.repl.lookupHelp(helpFunctionName, fullMode, arg);
            else
                contents = replab.infra.repl.callOriginalHelp(arg);
            end
        end
    else
        type = class(arg);
        if replab.compat.startsWith(type, 'replab')
            % RepLAB object
            contents = replab.infra.repl.lookupHelp(helpFunctionName, fullMode, type);
        else
            % Non RepLAB object, call original help function
            contents = replab.infra.repl.callOriginalHelp(arg);
        end
    end
    if replab.globals.runsInMatlabEmacs
        if ischar(arg)
            replab.infra.repl.matlabEmacsHelpDisplay(arg, contents);
        else
            replab.infra.repl.matlabEmacsHelpDisplay('object', contents);
        end
    else
        stdout = 1;
        fwrite(stdout, contents);
    end
end
