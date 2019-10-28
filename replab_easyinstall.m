function exitCode = replab_easyinstall
% Checks if the external dependencies are present, 
%
% If everything is setup correctly already, the script does not require user interaction, but will
% display output on the command line.
%
% If dependencies are missing, it will propose to the user to download and install the dependencies,
% from a ZIP file present in the main replab repository (in the ``externals`` branch).
% User input is required, in that case the script is interactive.
%
% This script asks user confirmation before doing anything.
%
%
% Returns:
%   integer: 0 if the external dependencies are present now, a positive number indicates an error
    
    isOctave = logical(exist('OCTAVE_VERSION', 'builtin'));
    scriptPath = mfilename('fullpath');
    [replabPath, ~, ~] = fileparts(scriptPath);
    
    externalPath = fullfile(replabPath, 'external');
    switch exist(externalPath)
      case 7
        % external directory exists, proceed
      case 0
        disp(['The external directory ' externalPath ' is not present.']);
        disp('This means that the VPI library is missing, and it is present by default in the RepLAB distribution.')
        disp('This script cannot help. Restart with a fresh copy of RepLAB');
        exitCode = 1;
      otherwise
        error('Path error for external'); 
    end

    
    MOcovPath = fullfile(externalPath, 'MOcov');
    MOcovPathExists = exist(MOcovPath) == 7;
    MOcovTestPath = fullfile(MOcovPath, 'README.md');
    MOcovZipPath = fullfile(externalPath, 'MOcov.zip');
    MOcovInstalled = MOcovPathExists && exist(MOcovTestPath) == 2;
    
    MOxUnitPath = fullfile(externalPath, 'MOxUnit');
    MOxUnitPathExists = exist(externalPath) == 7;
    MOxUnitTestPath = fullfile(MOxUnitPath, 'README.md');
    MOxUnitZipPath = fullfile(externalPath, 'MOxUnit.zip');
    MOxUnitInstalled = MOxUnitPathExists && exist(MOxUnitTestPath) == 2;
    
    SDPT3Path = fullfile(externalPath, 'SDPT3');
    SDPT3PathExists = exist(SDPT3Path) == 7;
    SDPT3TestPath = fullfile(SDPT3Path, 'README.md');
    SDPT3ZipPath = fullfile(externalPath, 'SDPT3.zip');
    SDPT3Installed = SDPT3PathExists && exist(SDPT3TestPath) == 2;
    
    YALMIPPath = fullfile(externalPath, 'YALMIP');
    YALMIPPathExists = exist(YALMIPPath) == 7;
    YALMIPTestPath = fullfile(YALMIPPath, 'README.txt');
    YALMIPZipPath = fullfile(externalPath, 'YALMIP.zip');
    YALMIPInstalled = YALMIPPathExists && exist(YALMIPTestPath) == 2;

    if MOcovInstalled && MOxUnitInstalled && SDPT3Installed && YALMIPInstalled
        disp('All external dependencies are present (or at least their README file).');
        disp('This script is thus not necessary');
        disp('You can directly run ''replab_addpaths'' to setup the RepLAB path.');
        exitCode = 0;
        return
    end
    
    % Test if we are in a Git repository
    gitPath = fullfile(replabPath, '.git');
    switch exist(gitPath)
      case 0
        % Does not contain a Git directory, proceed
      case 7
        % A Git directory exists, bail out
        disp('This script is designed to streamline the download of external dependencies.');
        disp(['However, this RepLAB directory ' replabPath ' is a Git repository, so the dependencies']);
        disp('should be updated using the Git submodule support.');
        disp('Please run ''git submodule init'' and ''git submodule update'' from the command line');
        disp('Or download and unzip the latest version without the Git data using the link below:');
        disp('https://github.com/replab/replab/archive/master.zip');
        disp('and run this script again.');
        exitCode = 1;
        return
      otherwise
        error([gitPath ' exists but is not a folder. Restart with a fresh copy of RepLAB.']);
    end
    
    disp('Some dependencies are missing. The script will now install the missing dependencies.');
    s = input('Do you wish to proceed [Y/n]?', 's');
    s = lower(strtrim(s));
    switch s
      case {'', 'y'}
        % proceed
      case 'n'
        exitCode = 1;
        return
      otherwise
        error('Unrecognized answer.');
    end
    
    zipName = 'replab_externals.zip';
    zipPath = fullfile(replabPath, zipName);
    url = 'https://github.com/replab/replab/raw/externals/replab_externals.zip';

    switch exist(zipPath)
      case 0
        disp('The ZIP file containing externals is not available locally.')
        disp('You have the option of downloading the file from');
        disp('https://github.com/replab/replab/raw/externals/replab_externals.zip');
        disp(['and place it in the root RepLAB folder ' replabPath]);
        disp('Or we can try to download the file automatically.');
        disp('After successful installation, this ZIP file can be deleted.');
        s = input('Do you wish to proceed with download [Y/n]?', 's');
        s = lower(strtrim(s));
        switch s
          case {'', 'y'}
            % proceed
          case 'n'
            exitCode = 1;
            return
          otherwise
            error('Unrecognized answer.');
        end
        try
            if isOctave
                [f, success] = urlwrite (url, zipPath);
                assert(success, 'Download did not work');
            else
                outfilename = websave(zipPath, url);
            end
        catch ME
            disp('Download failed. Download the file manually and restart the script.');
            exitCode = 1;
            return
        end            
        if exist(zipPath) ~= 2
            disp('Download failed, file not present. Download the file manually and restart the script.');
            exitCode = 1;
            return
        end
        s = dir(zipPath);
        if s.bytes < 1e6
            disp('Download failed, file is too small. Download the file manually and restart the script.');
            exitCode = 1;
            return
        end
      case 2
        % zip file exists, proceed
        disp(['Externals ZIP file ' zipPath ' available locally.']); 
        s = dir(zipPath);
        if s.bytes < 1e6
            disp('The external ZIP file is too small');
            exitCode = 1;
            return
        end
      otherwise
        error('Unrecognized type for the ZIP file path');
    end
    
    disp(['Unzipping ' zipPath]);
    try
        unzip(zipPath, externalPath);
    catch ME
        disp('Unzipping failed.')
        exitCode = 1;
        return
    end
    
    % MOcov
    if ~MOcovPathExists
        disp('Directory MOcovPath not present. Creating.');
        mkdir(externalPath, 'MOcov');
    end
    if MOcovInstalled
        disp('MOcov already present, skipping');
    else
        disp('Unzipping MOcov...');
        unzip(MOcovZipPath, MOcovPath);
    end
    
    % MOxUnit
    if ~MOxUnitPathExists
        disp('Directory MOxUnit not present. Creating.');
        mkdir(externalPath, 'MOxUnit');
    end
    if MOxUnitInstalled
        disp('MOxUnit already present, skipping');
    else
        disp('Unzipping MOxUnit...');
        unzip(MOxUnitZipPath, MOxUnitPath);
    end
    
    % SDPT3
    if ~SDPT3PathExists
        disp('Directory SDPT3 not present. Creating.');
        mkdir(externalPath, 'SDPT3');
    end
    if SDPT3Installed
        disp('SDPT3 already present, skipping');
    else
        disp('Unzipping SDPT3...');
        unzip(SDPT3ZipPath, SDPT3Path);
    end
    
    % YALMIP
    if ~YALMIPPathExists
        disp('Directory YALMIP not present. Creating.');
        mkdir(externalPath, 'YALMIP');
    end
    if YALMIPInstalled
        disp('YALMIP already present, skipping');
    else
        disp('Unzipping YALMIP...');
        unzip(YALMIPZipPath, YALMIPPath);
    end
    
    disp('');
    disp('Installation successful.');
    exitCode = 0;
end

