function exitCode = replab_easy_install
    isOctave = logical(exist('OCTAVE_VERSION', 'builtin'));
    scriptPath = mfilename('fullpath');
    [replabPath, ~, ~] = fileparts(scriptPath);
    externalPath = fullfile(replabPath, 'external');
    zipName = 'replab_externals.zip';
    zipPath = fullfile(replabPath, zipName);
    url = 'https://github.com/replab/replab/raw/externals/replab_externals.zip';
    MOcovPath = fullfile(externalPath, 'MOcov');
    MOcovTestPath = fullfile(MOcovPath, 'README.md');
    MOcovZipPath = fullfile(externalPath, 'MOcov.zip');
    MOxUnitPath = fullfile(externalPath, 'MOxUnit');
    MOxUnitTestPath = fullfile(MOxUnitPath, 'README.md');
    MOxUnitZipPath = fullfile(externalPath, 'MOxUnit.zip');
    SDPT3Path = fullfile(externalPath, 'SDPT3');
    SDPT3TestPath = fullfile(SDPT3Path, 'README.md');
    SDPT3ZipPath = fullfile(externalPath, 'SDPT3.zip');
    YALMIPPath = fullfile(externalPath, 'YALMIP');
    YALMIPTestPath = fullfile(YALMIPPath, 'README.md');
    YALMIPZipPath = fullfile(externalPath, 'YALMIP.zip');

    exitCode = 0;
    
    try
        if exist(zipPath) == 0
            disp('ZIP file containing externals unavailable, downloading...');
            if isOctave
                [f, success] = urlwrite (url, zipPath);
                assert(success, 'Download did not work');
            else
                outfilename = websave(zipPath, url);
            end
        else
            disp('ZIP file containing externals found');
        end
        assert(exist(zipPath) == 2);
        s = dir(zipPath);
        assert(s.bytes > 1e6, 'Downloaded file is too small');
        disp('Unzipping');
        unzip(zipPath, externalPath);
        
        % MOcov
        if exist(MOcovTestPath) == 2
            disp('MOcov already present, skipping');
        else
            disp('Unzipping MOcov...');
            unzip(MOcovZipPath, MOcovPath);
        end
        
        % MOxUnit

        if exist(MOxUnitTestPath) == 2
            disp('MOxUnit already present, skipping');
        else
            disp('Unzipping MOxUnit...');
            unzip(MOxUnitZipPath, MOxUnitPath);
        end
        
        % SDPT3

        if exist(SDPT3TestPath) == 2
            disp('SDPT3 already present, skipping');
        else
            disp('Unzipping SDPT3...');
            unzip(SDPT3ZipPath, SDPT3Path);
        end
        
        % YALMIP

        if exist(YALMIPTestPath) == 2
            disp('YALMIP already present, skipping');
        else
            disp('Unzipping YALMIP...');
            unzip(YALMIPZipPath, YALMIPPath);
        end
        
    catch ME
        ME
        ME.stack
        disp('Errored');
        exitCode = 1;
    end
    
    
    disp('Cleaning up ZIP files...')
    try
        delete(zipPath);
    catch ME
    end
    try
        delete(MOcovZipPath);
    catch ME
    end 
    try
        delete(MOxUnitZipPath);
    catch ME
    end
    try
        delete(SDPT3ZipPath);
    catch ME
    end
    try
        delete(YALMIPZipPath);
    catch ME
    end
end
