function help(varargin)
% Displays help on a function/object/...
%
% Args:
%   varargin: - if a single element, provides help for this element,
%               whether it is a variable, object, function, class,
%               package, ...
%             - if two elements, with one of them being the string '-f'
%               or '--full', provides full help on this element
%
% Note:
%   This function overloads matlab's 'help' function. Renaming this file
%   does not alter its functionality.
%
% Note:
%   The display of RepLAB objects is driven by the liquid-like templates
%   in ``src/+replab/+infra``. Those templates support Sphinx-like reference
%   syntax using backticks.

    % extract own name
    [pathStr, helpFunctionName, extension] = fileparts(which(mfilename));

    persistent codeBase

    if nargin == 1 && isequal(varargin{1}, '--clear')
        codeBase = [];
        return
    end

    % Are we in full mode or not?
    fullMode = false;
    if (length(varargin) == 2)
        if isequal(varargin{1}, '-f') || isequal(varargin{1}, '--full')
            fullMode = true;
            varargin = {varargin{2}};
        elseif isequal(varargin{2}, '-f') || isequal(varargin{2}, '--full')
            fullMode = true;
            varargin = {varargin{1}};
        end
    end

    % Are we asking help on an object or method applied to an object?
    if (length(varargin) == 1) && ischar(varargin{1})
        % we check if the subject starts with a variable
        dotPositions = strfind(varargin{1}, '.');
        if isempty(dotPositions)
            firstPart = varargin{1};
            secondPart = '';
        else
            firstPart = varargin{1}(1:dotPositions(1)-1);
            secondPart = varargin{1}(dotPositions(1)+1:end);
        end
        isVar = 0;
        try
            isVar = evalin('caller', ['exist(''', firstPart, ''')']);
        catch
        end
        if isVar == 1
            % object of interest is a variable or object (or method applied
            % to an object)
            type = evalin('caller', ['class(', firstPart, ')']);
            if (length(type) >= 6) && isequal(type(1:6), 'replab')
                % object is a replab object
                if isempty(secondPart)
                    varargin{1} = type;
                    msg = sprintf('%s is an object of type %s.', firstPart, ...
                                  replab.infra.repl.linkHelp(helpFunctionName, type, type));
                    disp(msg);
                    disp(' ');
                else
                    varargin{1} = [type, '.', secondPart];
                end
                msg = sprintf('--- help for %s ---', replab.infra.strong(varargin{1}));
                disp(msg);
            else
                % non-replab variable/object. We need to make a copy to our
                % workspace for the next function to see it
                eval([firstPart, ' = evalin(''caller'', ''', firstPart, ''');']);
            end
        end
    end

    if (length(varargin) == 1) && (length(varargin{1}) >= 6) && (isequal(varargin{1}(1:6), 'replab')) ...
            && ((length(varargin{1}) < 7) || (varargin{1}(7) == '.'))
        % We are looking for a replab-related help
        name = varargin{1};
        if isempty(codeBase)
            fprintf('Building help index...');
            [srcRoot, ~, ~] = fileparts(mfilename('fullpath'));
            codeBase = replab.infra.crawl(srcRoot);
            disp('done.');
        end
        parts = strsplit(name, '.');
        el = codeBase.get(parts{:});
        hasToggle = false;
        switch class(el)
          case 'replab.infra.Package'
            hasToggle = false;
            templateSuffix = 'package';
            docEl = el;
          case 'replab.infra.Class'
            hasToggle = true;
            templateSuffix = 'class';
            docEl = el;
          case 'replab.infra.Function'
            hasToggle = true;
            templateSuffix = 'function';
            docEl = el;
          case {'replab.infra.ConcreteClassElement', 'replab.infra.InheritedClassElement'}
            hasToggle = true;
            templateSuffix = 'class_element';
            docEl = el.declarations.findBestDocumented;
          otherwise
            error('replab:helpError', 'Object referenced by %s is of type %s', name, class(element));
        end
        if ~isempty(docEl) && (~isa(docEl, 'replab.infra.SourceElement') || docEl.doc.isempty)
            docEl = [];
        end
        if fullMode
            helpCurrent = [helpFunctionName ' -f'];
            helpToggled = helpFunctionName;
        else
            helpCurrent = helpFunctionName;
            helpToggled = [helpFunctionName ' -f'];
        end
        if hasToggle && fullMode
            templateName = ['help_full_' templateSuffix];
        else
            templateName = ['help_' templateSuffix];
        end
        stdout = 1;
        fullId = el.fullIdentifier;
        fwrite(stdout, replab.infra.repl.templateHelp(templateName, el, docEl, helpCurrent, {fullId}, {fullId}));
        if hasToggle && replab.settings.consoleUseHTML
            if fullMode
                link = replab.infra.repl.linkHelp(helpToggled, 'help mode', fullId);
                disp(['   Toggle display to ', link]);
            else
                link = replab.infra.repl.linkHelp(helpToggled, 'reference mode', fullId);
                disp(['   Toggle display to ', link]);
            end
        end
        if isa(el, 'replab.infra.SourceElement')
            disp(' ');
            disp(['   ', replab.infra.repl.linkOpen('See source', '', el.absoluteFilename, el.startLineNumber)]);
        end
        fprintf('\n\n');
%        separationPattern = '   .  ';
%        disp(repmat(separationPattern, 1, floor(replab.settings.strMaxColumns/length(separationPattern))));
%        disp(' ');
    else
        if isempty(replab.settings.systemHelpPath)
            try
                % Some clear happened, we try to capture matlab's help
                % again...
                replab_init(0);
            catch
                error('The matlab help path was not captured. Please use replab_init first.');
            end
        end

        % We call matlab's help function
        currentPath = strrep(pwd, '\', '/');

        if ~replab.compat.isOctave
            cd(replab.settings.systemHelpPath);
            message = [];
            try
                help(varargin{:});
            catch message
            end
            cd(currentPath);
            if ~isempty(message)
                error(message);
            end
        else
            % In some versions of octave earlier than 5.1.0, the
            % current path had a lower priority than the path order.
            % Then we also need replab's path...

            replabHelpPath = fileparts(which('replab_init'));
            replabHelpPath = [strrep(replabHelpPath, '\', '/'), '/src'];

            cd(replab.settings.systemHelpPath);
            addpath(replab.settings.systemHelpPath);
            message = [];
            try
                help(varargin{:});
            catch message
            end
            cd(currentPath);
            addpath(replabHelpPath);
            if ~isempty(message)
                error(message);
            end
        end
    end
end
