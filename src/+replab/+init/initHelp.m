function initHelp(verbose)
% Initializes the RepLAB help system
%
% In particular, it attemps to locate the original Matlab/Octave help function and preserve a handle of it;
% if the help function is not already shadowed, it shadows it with the integrated RepLAB help system
%
% Args:
%   verbose ({0, 1, 2}): Controls the display level

    candidates = replab.init.findInstancesInPath('help'); % Finds all candidates in the path for "help"
    basePath = replab.globals.replabPath;
    enableHelpOverload = true;
    matlabEmacsHelpPath = [];
    originalHelpPath = [];
    replabHelpPath = [];
    unknownHelpPaths = {};
    targetPatternMatlab = 'Copyright 1984-\d\d\d\d The MathWorks, Inc.';
    targetPatternOctave = 'This file is part of Octave.';
    targetPatternEmacs = 'EMACSCAP';
    for i = 1:length(candidates)
        candidate = strrep(candidates{i}, '\', '/');
        if (replab.compat.endsWith(candidate, 'toolbox/matlab/helptools/help.m') ...
                || ~isempty(regexp(fileread(candidate), targetPatternMatlab))) ...
                && isempty(originalHelpPath)
            if replab.compat.isOctave
                warning('It looks like Matlab''s help function is in the path while using Octave');
                enableHelpOverload = false;
            end
            originalHelpPath = candidate;
        elseif (replab.compat.endsWith(candidate, 'm/help/help.m') ...
                || ~isempty(regexp(fileread(candidate), targetPatternOctave))) ...
                && isempty(originalHelpPath)
            if ~replab.compat.isOctave
                warning('It looks like Octave''s help function is in the path while using Matlab');
                enableHelpOverload = false; % same
            end
            originalHelpPath = candidate;
        elseif replab.compat.endsWith(candidate, 'toolbox/help.m') ...
                && ~isempty(strfind(fileread(candidate), targetPatternEmacs)) ...
                && isempty(matlabEmacsHelpPath)
            matlabEmacsHelpPath = candidate;
        elseif replab.compat.endsWith(candidate, '+replab/src/help_overload/help.m') && isempty(replabHelpPath)
            replabHelpPath = candidate;
        else
            unknownHelpPaths{end+1} = candidate;
        end
    end
    myHelpPath = fullfile(replab.globals.replabPath, 'src', 'help_overload', 'help.m');

    if ~isempty(replabHelpPath) && ~isequal(replabHelpPath, myHelpPath)
        error('Another version %s of the RepLAB help overload %s already exists in the path.', ...
              replabHelpPath, myHelpPath)
    end
    if ~isempty(unknownHelpPaths)
        warning('Several unidentified help functions present in path: %s', strjoin(unknownHelpPaths, pathsep))
        enableHelpOverload = false;
    end
    if enableHelpOverload
        if isempty(originalHelpPath)
            warning('Matlab/Octave original help function not present in path');
            handle = @(arg) 'Help unavailable'
        else
            currentPathStr = pwd;
            try
                cd(originalHelpPath(1:end-6));
                handle = @help;
            catch
                warning('Error when capturing original help function');
                enableHelpOverload = false;
            end
            cd(currentPathStr);
        end
        replab.globals.defaultHelpFunction(handle);
        addpath(fullfile(basePath, 'src', 'help_overload'));
    else
        disp('To enable help functionality, copy the file src/help_overload/help.m to src/replab_help.m');
    end
    if ~isempty(matlabEmacsHelpPath)
        replab.globals.runsInMatlabEmacs(true);
    end
end
