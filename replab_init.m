function replab_init(varargin)
% function replab_init
%
% Sets up the search path in order to enable all functionalities of the RepLAB library.
%
% Also verifies that a SDP solver is available, installs and registers the bundled SDPT3 solver if needed,
% and initializes important variables.
%
% The global variable `+replab.+globals.replabPath` is set up and ``mlock``-ed so that it persists a ``clear all``.
%
% This script accepts the following options:
%
% * The key ``verbose`` followed by a value 0,1,2, controls the display level when initializating
%
%   * 0: only produce a display in case of error/warning or for critical cases
%   * 1: informs of the changes made (default value)
%   * 2: prints full information
%
% * The option ``autoinstall``: Install dependencies automatically
%
% * The option ``showurls``: Shows the URL that need to be downloaded to activate dependencies
%
% Args:
%     varargin (cell(1,\*) of key/value pairs): Options, see description above
%
% Example:
%     >>> replab_init   % doctest: +SKIP
%         Adding RepLAB to the path
%         Adding VPI to the path
%         Adding MOxUnit to the path
%         Adding embedded YALMIP to the path
%         Adding MOcov to the path

    persistent allGood % set to true if everything is already set up previously, including running ``replab_config``

    verbose = 1;
    autoinstall = false;
    showurls = false;
    i = 1;
    while i <= nargin
        key = varargin{i};
        i = i + 1;
        switch key
          case 'showurls'
            showurls = true;
          case 'autoinstall'
            autoinstall = true;
          case 'verbose'
            value = varargin{i};
            if isa(value, 'char')
                value = str2num(value);
            end
            assert(isa(value, 'double'));
            verbose = value;
            i = i + 1;
          otherwise
            error('Unrecognized option %s', key);
        end
    end


    %% If everything was already set up previously we exit rapidly
    if isempty(allGood)
        allGood = false;
    end
    if allGood && ~showurls
        replab.init.log_(2, 'replab_init has already been successfully called.');
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

    % Check if the current instance of RepLAB is already in the path, and
    % throw an error if another instance of RepLAB is already in the path)
    alreadyInPath = isInPath('replab_init', '', basePath, verbose);
    if ~alreadyInPath
        addpath(basePath);
        if verbose >= 1
            disp('Adding RepLAB to the path');
        end
    end

    srcAlreadyInPath = isInPath('replab_Version', 'src', basePath, verbose);
    if ~srcAlreadyInPath
        addpath(fullfile(basePath, 'src'));
        if (verbose == 1) && (alreadyInPath)
            % At verbose level 1 we do not make a difference between the
            % replab folder and subfolder, adding any one of them leads
            disp('Adding RepLAB to the path');
        elseif verbose >= 2
            disp('Adding RepLAB src folder to the path');
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

    %% Memorizing the options
    replab.globals.verboseInit(verbose);
    replab.globals.autoInstall(autoinstall);

    %% Run the config script
    if showurls
        replab.init.showUrls;
    else
        replab_config
        % If we haven't errored yet, then it's all good
        allGood = true;
    end

end


function alreadyInPath = isInPath(functionName, subfolder, basePath, verbose)
% Checks the presence of a RepLAB function in the path
%
% Throws an error if a function with the resired name is found in another
% place than the current RepLAB library folder. This guarantees that no
% other instance of RepLAB is also present in the path.
%
% Args:
%   functionName (charstring): Name of the RepLAB function we are looking
%                              for
%   subfolder (charstring): RepLAB subfolder in which the function should
%                           be found
%   basePath (charstring): base path of the current RepLAB library
%   isOctave (logical): true if the platform is octave
%   verbose ({0, 1, 2}): Controls the display level
%
% Returns:
%   logical: whether the path already contains the

    allPaths = strsplit(path, pathsep);
    candidates = findInstancesInPath(functionName);
    versionFilePath = strrep(fullfile(basePath, subfolder, [functionName, '.m']), '\', '/');
    alreadyInPath = false;
    for i = 1:length(candidates)
        candidate = strrep(candidates{i}, '\', '/');
        if isequal(candidate, versionFilePath)
            alreadyInPath = true;
            if verbose >= 2
                if isempty(subfolder)
                    disp('RepLAB is already in the path');
                else
                    disp('RepLAB subfolder is already in the path');
                end
            end
        else
            error(['Another instance of RepLAB in folder ' fileparts(fileparts(candidate)) ' ' ...
                   'is already in the path' char(10) ...
                   'Use this one or remove it from the path.']);
        end
    end
end

function str = findInstancesInPath(item)
% Finds all implementations of a function/class in the current path
%
% Returns the same results as ``which(item, '-all')``, except it only considers
% the current path if it is explicitly present .
%
% Assumes that ``item`` is implemented as a ``.m`` file.
%
% Note: this is a copy of `+replab.+init.findInstancesInPath`
%
% Args:
%   item (charstring): MATLAB function to look for
%
% Returns:
%   cell(\*,1) of charstring: Paths to the code files
    str = {};
    paths = strsplit(path, pathsep);
    for i = 1:length(paths)
        if ~isequal(paths{i}, '.') % Octave has '.' in the path
            candidate = fullfile(paths{i}, [item '.m']);
            if exist(candidate) == 2
                str{end+1,1} = candidate;
            end
        end
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
