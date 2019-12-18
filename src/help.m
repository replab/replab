function help(varargin)
% Displays help on a function/object/...
%
% Arg:
%     varargin: - if a single element, provides help for this element,
%                 whether it is a variable, object, function, class, 
%                 package, ...
%               - if two elements, with one of them being the string '-f'
%                 or '--full', provides full help on this element
%
% Note:
%     This function overloads matlab's 'help' function. Renaming this file
%     does not alter its functionality.

    % extract own name
    [pathStr, helpFunctionName, extension] = fileparts(which(mfilename));

    persistent codeBase
    
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
                    replab.infra.dispH([firstPart, ' is an object of type ', type, '.'], firstPart, helpFunctionName, fullMode);
                    disp(' ');
                else
                    varargin{1} = [type, '.', secondPart];
                end
                replab.infra.dispH(['--- help for ', varargin{1}, ' ---'], varargin{1}, helpFunctionName, fullMode);
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
            codeBase = replab.infra.CodeBase.crawl(srcRoot);
            disp('done.');
        end
        parts = strsplit(name, '.');
        element = codeBase.get(parts{:});
        switch class(element)
          case 'replab.infra.Package'
            help_package(codeBase, element, helpFunctionName, fullMode);
          case 'replab.infra.Class'
            help_class(codeBase, element, helpFunctionName, fullMode);
          case 'replab.infra.Function'
            help_function(codeBase, element, helpFunctionName, fullMode);
          case {'replab.infra.ConcreteClassElement', 'replab.infra.InheritedClassElement'}
            help_classElement(codeBase, element, helpFunctionName, fullMode);
          otherwise
            error('replab:helpError', 'Object referenced by %s is of type %s', name, class(element));
        end
    else
        if isempty(replab.Parameters.matlabHelpPath)
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

        if ~replab.platformIsOctave
            cd(replab.Parameters.matlabHelpPath);
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

            cd(replab.Parameters.matlabHelpPath);
            addpath(replab.Parameters.matlabHelpPath);
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

function help_package(codeBase, package, helpFunctionName, fullMode)
    fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_package.liquid');
    tmpl = replab.lobster.Template.load(fn);
    disp(tmpl.render(struct('package', package)));
end


function help_class(codeBase, element, helpFunctionName, fullMode)
    if fullMode
        fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_class.liquid');
    else
        fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_full_class.liquid');
    end
    tmpl = replab.lobster.Template.load(fn);
    disp(tmpl.render(struct('cls', element)));
    if replab.Parameters.consoleUseHTML
        if fullMode
            link = replab.infra.linkHelp(helpFunctionName, element.fullIdentifier, element.fullIdentifier);
            disp(sprintf('Toggle display to help page for %s', link));
        else
            link = replab.infra.linkHelp(helpFunctionName, element.fullIdentifier, element.fullIdentifier, '-f');
            disp(sprintf('Toggle display to reference page for %s', link));
        end
        disp(replab.infra.linkOpen('See source', '', element.absoluteFilename, element.startLineNumber));
    end
end

function help_classElement(codeBase, element, helpFunctionName, fullMode)
    if fullMode
        fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_class_element.liquid');
    else
        fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_full_class_element.liquid');
    end
    tmpl = replab.lobster.Template.load(fn);
    disp(tmpl.render(struct('el', element)));
    if replab.Parameters.consoleUseHTML
        if fullMode
            link = replab.infra.linkHelp(helpFunctionName, element.fullIdentifier, element.fullIdentifier);
            disp(sprintf('Toggle display to help page for %s', link));
        else
            link = replab.infra.linkHelp(helpFunctionName, element.fullIdentifier, element.fullIdentifier, '-f');
            disp(sprintf('Toggle display to reference page for %s', link));
        end
        disp(replab.infra.linkOpen('See source', '', element.absoluteFilename, element.startLineNumber));
    end
end

function help_function(codeBase, fun, helpFunctionName, fullMode)
    if fullMode
        fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_function.liquid');
    else
        fn = fullfile(codeBase.rootFolder, '+replab', '+infra', 'help_full_function.liquid');
    end
    tmpl = replab.lobster.Template.load(fn);
    disp(tmpl.render(struct('fun', fun)));
    if replab.Parameters.consoleUseHTML
        if fullMode
            link = replab.infra.linkHelp(helpFunctionName, fun.fullIdentifier, fun.fullIdentifier);
            disp(sprintf('Toggle display to help page for %s', link));
        else
            link = replab.infra.linkHelp(helpFunctionName, fun.fullIdentifier, fun.fullIdentifier, '-f');
            disp(sprintf('Toggle display to reference page for %s', link));
        end
        disp(replab.infra.linkOpen('See source', '', fun.absoluteFilename, fun.startLineNumber));
    end
end
