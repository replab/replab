function replab_addpaths(verbose)
    % replab_addpaths([verbose])
    %
    % replab_addpaths adds directories that are needed to the search path
    % in order to enable all functionalities of the library.
    %
    % The optional parameter 'verbose' controls the display level:
    %   0 only produce a display in case of error
    %   1 informs of the changes made (default value)
    %   2 prints the complete information
    
    
    %% Parameter checking
    if nargin < 1
        verbose = 1;
    end
    
    
    %% Action
    [pathStr, name, extension] = fileparts(which(mfilename));
    addpath(pathStr);
    
    % Making sure the VPI library is in the path and working
    VPIInPath = false;
    try
        VPIInPath = isequal(num2str(vpi('1234567890098765432100123456789')), '    1234567890098765432100123456789');
    catch
    end
    if ~VPIInPath
        if exist([pathStr '/external/vpi/@vpi/vpi.m']) ~= 2
            error(['The VPI library was not found in the folder ', pathStr, '/external/vpi']);
        else
            addpath([pathStr '/external/vpi']);
            if verbose >= 1
                disp('Adding VPI to the path');
            end
        end
    elseif verbose >= 2
        disp('VPI is already in the path');
    end
    
    % Making sure MOxUnit is in the path and working
    MOxUnitInPath = false;
    try
        moxunit_set_path;
        MOxUnitInPath = true;
    catch
    end
    if ~MOxUnitInPath
        if exist([pathStr '/external/MOxUnit/MOxUnit/moxunit_set_path.m']) ~= 2
            warning(['The MOxUnit library was not found in the folder ', pathStr, '/external/MOxUnit', char(10), ...
                'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                'It will not be possible to run the tests.']);
        else
            addpath([pathStr '/external/MOxUnit/MOxUnit']);
            moxunit_set_path;
            if verbose >= 1
                disp('Adding MOxUnit to the path');
            end
        end
    elseif verbose >= 2
        disp('MOxUnit is already in the path');
    end
    
    % Making sure MOcov is in the path and working
    MOcovInPath = false;
    try
        mocov_get_absolute_path('.');
        MOcovInPath = true;
    catch
    end
    if ~MOcovInPath
        if exist([pathStr '/external/MOcov/MOcov/mocov.m']) ~= 2
            warning(['The MOcov library was not found in the folder ', pathStr, '/external/MOcov', char(10), ...
                'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                'It will not be possible to run the tests.']);
        else
            addpath([pathStr '/external/MOcov/MOcov']);
            moxunit_set_path;
            if verbose >= 1
                disp('Adding MOcov to the path');
            end
        end
    elseif verbose >= 2
        disp('MOcov is already in the path');
    end

    
    % Making sure YALMIP is in the path and working
    YALMIPInPath = false;
    try
        yalmip('version');
        YALMIPInPath = true;
    catch
    end
    if ~YALMIPInPath
        if verbose >= 1
            warning('YALMIP was not found in the path, some functionalities of the library might be disabled');
        end
    elseif verbose >= 2
        disp('YALMIP is already in the path');
    end
end
