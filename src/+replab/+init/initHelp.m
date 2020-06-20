function res = initHelp(verbose)
% Initializes the RepLAB help system
%
% In particular, it attemps to locate the original Matlab/Octave help function and preserve a handle of it;
% if the help function is not already shadowed, it shadows it with the integrated RepLAB help system

    % Finds all candidates in the path for "help"
    if isOctave
        helpCandidates = {};
        paths = strsplit(path, pathsep);
        for i = 1:length(paths)
            candidate = fullfile(paths{i}, 'help.m');
            if exist(candidate) == 2
                helpCandidates{end+1} = candidate;
            end
        end
    else
        helpCandidates = which('help', '-all');
    end


        if enableHelpOverload
        replab.globals.defaultHelpFunction(handle);
        addpath(fullfile(pathStr, 'src', 'help'));
    end

    if ~isempty(matlabEmacsHelpPath)
        replab.globals.runsInMatlabEmacs(true);
    end


    %% Action -- first capture the folder containg matlab's help.m

    enableHelpOverload = true;
    if 1 % isempty(replab.globals.defaultHelpFunction) resolve this later
        matlabEmacsHelpPath = [];
        originalHelpPath = [];
        replabHelpPaths = {};
        unknownHelpPaths = {};
        for i = 1:length(helpCandidates)
            candidate = helpCandidates{i};
            if ~isempty(regexp(candidate, 'toolbox[/\\]matlab[/\\]helptools[/\\]help\.m$')) ...
                    && isempty(originalHelpPath)
                if isOctave
                    warning('It looks like Matlab''s help function is in the path while using Octave');
                    enableHelpOverload = false; % we shouldn't have Matlab help in path when using Octave
                end
                originalHelpPath = candidate;
            elseif ~isempty(regexp(candidate, 'm[/\\]help[/\\]help\.m$')) ...
                    && isempty(originalHelpPath)
                if ~isOctave
                    warning('It looks like Octave''s help function is in hte path while using Matlab');
                    enableHelpOverload = false; % same
                end
                originalHelpPath = candidate;
            elseif ~isempty(regexp(candidate, 'toolbox[\\/]help\.m$')) ...
                    && ~isempty(strfind(fileread(candidate), 'EMACSCAP'))  && isempty(matlabEmacsHelpPath)
                matlabEmacsHelpPath = candidate;
            elseif ~isempty(regexp(candidate, 'src[\\/]help[\\/]help\.m$')) ...
                    && ~isempty(strfind(fileread(candidate), 'replab'))
                replabHelpPaths{end+1} = candidate;
            else
                unknownHelpPaths{end+1} = candidate;
            end
        end
        if ~isempty(replabHelpPaths)
            error('A version of help already exists in the path.')
            % TODO: handle when I'm myself
        end
        if ~isempty(unknownHelpPaths)
            warning(['Several unidentified help functions present in path: ' strjoin(unknownHelpPaths, pathsep)]);
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
        end
    end
end
