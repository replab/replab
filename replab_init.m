function replab_init(verbose)
% Sets up the search path in order to enable all functionalities of the RepLAB library.
%
% Also verifies that a SDP solver is available, installs and registers the bundled SDPT3 solver if needed,
% and initializes important variables.
%
% A few persistent variables in functions present in `replab.globals` are initialized by this script;
% those persistent variables are locked by `mlock` to avoid being cleared when ``clear all`` is used;
% this concerns in particular `replab.globals.replabPath`, `replab.globals.defaultHelpFunction` and
% `replab.globals.runsInMatlabEmacs`.
%
% Args:
%     verbose ({0, 1, 2}, optional): Controls the display level
%                                    - 0: only produce a display in case of error/warning or for critical cases
%                                    - 1: informs of the changes made (default value)
%                                    - 2: prints full information
%
% Example:
%     >>> replab_init   % doctest: +SKIP
%         Adding RepLAB to the path
%         Adding RepLAB package to the path
%         Adding VPI to the path
%         Adding MOxUnit to the path
%         Adding embedded YALMIP to the path
%         Adding MOcov to the path

    persistent allGood % set to true if everything is already set up previously

    %% Parameter checking
    if nargin < 1
        verbose = 1;
    end


    %% If everything was already set up previously we exit rapidly
    if isempty(allGood)
        allGood = false;
    end
    if allGood
        if verbose >= 2
            disp('Replab_init has already been successfully called.');
        end
        return
    end

    %% Let us check the Matlab/Octave version
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if ~isOctave
        platform = 'Matlab';
        minimalVersion = '8.6';
        currentVersion = version;
        currentVersion = currentVersion(1:find(currentVersion==' ',1)-1);
    else
        platform = 'Octave';
        minimalVersion = '4.2.2';
        currentVersion = version;
    end
    if ~isLaterVersion(minimalVersion, currentVersion)
        warning(['Current version of ', platform, ' is ', currentVersion, ' but the minimal supported version is ', minimalVersion, '.']);
    end

    %% Adding the RepLAB source folder to the path
    [basePath, name, extension] = fileparts(mfilename('fullpath'));
    basePath = strrep(basePath, '\', '/'); % we normalize to Unix style filesep, as it is compatible with all

    % Check if another instance of RepLAB is already in the path
    candidates = whichAll('replab_Version');
    versionFilePath = strrep(fullfile(basePath, 'src', 'replab_Version.m'), '\', '/');
    alreadyInPath = false;
    for i = 1:length(candidates)
        candidate = strrep(candidates{i}, '\', '/');
        if isequal(candidate, versionFilePath)
            alreadyInPath = true;
            if verbose >= 2
                disp('RepLAB is already in the path');
            end
        else
            error(['Another instance of RepLAB in folder ' fileparts(fileparts(candidate)) ...
                   'is already in the path' char(10) 'Use this one or remove it from the path.']);
        end
    end

    if ~alreadyInPath
        addpath(fullfile(basePath, 'src'));
        if verbose >= 1
            disp('Adding RepLAB to the path');
        end
    end


    %% Memorizing RepLAB root folder if not done before
    rgrp = replab.globals.replabPath;
    if isempty(rgrp)
        replab.globals.replabPath(basePath);
    elseif ~isequal(rgrp, basePath)
        error('Mismatch between memorized RepLAB path %s and current RepLAB path %s', rgrp, basePath);
    else
        % path already memorized
    end

    % From now on, we can use stuff in the replab.init package

    %% Initializes the RepLAB help system
    replab.init.initHelp(verbose);

    %% Verifying that the nlinfit function is available
    replab.init.initNlinfit(verbose);

    %% Verifying that the symbolic toolbox is available
    replab.init.initSym(verbose);

    %% VPI
    VPIInPath = replab.init.initVPI(verbose);

    %% MOxUnit
    MOxUnitInPath = replab.init.initMOxUnit(verbose);

    %% YALMIP
    YALMIPInPath = replab.init.initYALMIP(verbose);

    %% SDPT3
    % Makes sure a working SDP solver is in the path and working, otherwise tries to add SDPT3
    SDPSolverInPath = false;
    if YALMIPInPath
        SDPSolverInPath = replab.init.initSDP(verbose);
    elseif verbose >= 2
        warning('YALMIP not in path, not checking for an SDP solver')
    end

    %% MOcov
    MOcovInPath = replab.init.initMOcov(verbose);

    %% If everything was successful, the next call will be quicker
    if VPIInPath && MOxUnitInPath && YALMIPInPath && SDPSolverInPath && MOcovInPath
        allGood = true;
    end

    function str = whichAll(item)
    % Implementation of the which function for Octave compatibility
    %
    % Returns the same results as ``which(item, '-all')``.
    %
    % Assumes that ``item`` is implemented as a ``.m`` file.
    %
    % Note: this is a copy of replab.compat.whichAll
    %
    % Args:
    %   item (charstring): MATLAB function to look for
    %
    % Returns:
    %   cell(\*,1) of charstring: Paths to the code files
        if isOctave
            str = {};
            paths = strsplit(path, pathsep);
            for i = 1:length(paths)
                candidate = fullfile(paths{i}, [item '.m']);
                if exist(candidate) == 2
                    str{end+1,1} = candidate;
                end
            end
        else
            str = which(item, '-all');
        end
    end

    function ok = isLaterVersion(minimalVersion, currentVersion)
    % isLaterVersion - checks if the current version is old enough
    %
    % ok = isLaterVersion(minimalVersion, currentVersion)
    %
    % Args:
    %     minimalVersion (charstring) : threshold version
    %     currentVersion (charstring) : actual version
    %
    % Returns:
    %     logical: answer
    %
    % Example:
    %     isLaterVersion('4.2.2', '04.3') % true
    %     isLaterVersion('4.2.2', '4.2.2.01') % true

        minimalVersion = [minimalVersion, '.'];
        currentVersion = [currentVersion, '.'];

        while (sum(minimalVersion == '.') >= 1) && (sum(currentVersion == '.') >= 1)
            minPoint = find(minimalVersion == '.', 1);
            currPoint = find(currentVersion == '.', 1);
            if str2num(minimalVersion(1:minPoint-1)) < str2num(currentVersion(1:currPoint-1))
                ok = true;
                return;
            elseif str2num(minimalVersion(1:minPoint-1)) > str2num(currentVersion(1:currPoint-1))
                ok = false;
                return;
            end
            minimalVersion = minimalVersion(minPoint+1:end);
            currentVersion = currentVersion(currPoint+1:end);
        end

        ok = ~((sum(minimalVersion == '.') >= 1) && (sum(currentVersion == '.') < 1));
    end

end
