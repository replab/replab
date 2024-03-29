function [contents, docTopic] = help(varargin)
% Displays help on a function/object/...
%
% Args:
%   varargin: - if a single element, provides help for this element,
%               whether it is a variable, object, function, class,
%               package, ...
%             - if two elements, with one of them being the string '-f'
%               or '-full', provides full help on this element
%
% Returns
% -------
%   contents: charstring
%     The help information in text format.
%   docTopic: charstring
%     Name of the associated reference page (for matlab internal functions
%     only)
%
% Note:
%   This function overloads matlab's 'help' function. Renaming this file
%   does not alter its functionality.
%
% Note:
%   The display of RepLAB objects is driven by the liquid-like templates
%   in ``src/+replab/+infra``. Those templates support Sphinx-like reference
%   syntax using backticks.

    contents = '';
    docTopic = '';

    % This function is part of replab, so it should only be called once the
    % replab package has been initialized. We provide a warning message if
    % this does not seem to be the case.
    repPath = '';
    try
        repPath = replab.globals.replabPath;
    catch
    end
    if isempty(repPath)
        pathToHelpFolder = fileparts(mfilename('fullpath'));        
        warning(['A RepLAB function is being called, but it seems that the RepLAB library was not initialized. Make sure that the folder ', pathToHelpFolder, ' was not inadvertantly included in your matlab path.']);
    end
    try
        replab_init;
    catch
        pathToHelpFolder = fileparts(mfilename('fullpath'));
        replabPos = strfind(lower(pathToHelpFolder), 'replab');
        if isempty(replabPos)
            pathToReplab = 'replab';
        else
            pathToReplab = pathToHelpFolder(1:replabPos+5);
        end
        error(['A RepLAB function is being called, but it seems that the RepLAB library cannot be initialized. Add the ', pathToReplab, ' folder in your path if you want to use this library, or make sure that the folder ', pathToHelpFolder, ' was not inadvertantly included in your matlab path.']);
    end
    
    if nargin == 1 && isequal(varargin{1}, '--clear')
        replab.globals.codeBase([]);
        return
    end

    % Are we in full mode or not?
    fullMode = false;
    if nargin >= 1
        arg = varargin{1};
    else
        arg = [];
    end
    if nargin == 2
        if isequal(varargin{1}, '-f') || isequal(varargin{1}, '--full')
            fullMode = true;
            arg = varargin{2};
        elseif isequal(varargin{2}, '-f') || isequal(varargin{2}, '--full')
            fullMode = true;
        end
    end

    % extract own name
    [~, helpFunctionName, ~] = fileparts(which(mfilename));
    helpPrinted = false;
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
                helpPrinted = (nargout==0) && (~replab.globals.runsInMatlabEmacs) && (~replab.compat.isOctave);
                [contents, docTopic] = replab.infra.repl.callOriginalHelp({varargin}, helpPrinted, varName, varValue);
            end
        else
            if replab.compat.startsWith(arg, 'replab') && ~isequal(arg, 'replab_init') && ~isequal(arg, 'replab_easyinstall')
                % it is a RepLAB command
                contents = replab.infra.repl.lookupHelp(helpFunctionName, fullMode, arg);
            else
                helpPrinted = (nargout==0) && (~replab.globals.runsInMatlabEmacs) && (~replab.compat.isOctave);
                [contents, docTopic] = replab.infra.repl.callOriginalHelp({varargin}, helpPrinted);
            end
        end
    else
        type = class(arg);
        if replab.compat.startsWith(type, 'replab')
            % RepLAB object
            contents = replab.infra.repl.lookupHelp(helpFunctionName, fullMode, type);
        else
            % Non RepLAB object, call original help function
            helpPrinted = (nargout==0) && (~replab.globals.runsInMatlabEmacs) && (~replab.compat.isOctave);
            [contents, docTopic] = replab.infra.repl.callOriginalHelp({varargin}, helpPrinted);
        end
    end
    if nargout == 0
        if ~helpPrinted
            % The help was not displayed yet, we do it now
            if replab.globals.runsInMatlabEmacs
                if ischar(arg)
                    replab.infra.repl.matlabEmacsHelpDisplay(arg, contents);
                else
                    replab.infra.repl.matlabEmacsHelpDisplay('object', contents);
                end
            else
                stdout = 1;
                fwrite(stdout, [char(10), contents, char(10)]);
            end
        end
        clear contents; % to avoid displaying twice if called from the command line
    end
end
