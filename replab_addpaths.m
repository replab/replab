function replab_addpaths(verbose, includeEmbeddedSolver)
    % replab_addpaths([verbose], [includeEmbeddedSolver])
    %
    % replab_addpaths adds directories that are needed to the search path
    % in order to enable all functionalities of the library.
    %
    % The optional parameter 'verbose' controls the display level:
    %   0 only produce a display in case of error/warning or for critical
    %     cases
    %   1 informs of the changes made (default value)
    %   2 prints the complete information
    %
    % The optional parameter 'includeEmbeddedSolver' is 0 by default. Set
    % it to 1 to add the folders 'external/YALMIP/...' and 'externl/SDPT3'
    % to the path (avoid this if you already have your own installation of
    % these libraries)
    %
    % example: replab_addpaths
    
    %% Parameter checking
    if nargin < 1
        verbose = 1;
    end

    if nargin < 2
        includeEmbeddedSolver = 0;
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
    
    if includeEmbeddedSolver == 1
        % Making sure YALMIP is in the path and working
        YALMIPInPath = false;
        try
            yalmip('version');
            YALMIPInPath = true;
        catch
        end
        if ~YALMIPInPath
            if exist([pathStr '/external/YALMIP/yalmiptest.m']) ~= 2
                warning(['The YALMIP library was found neither in the path nor in the folder ', pathStr, '/external/YALMIP, some functionalities of the library might be disabled']);
            else
                addpath([pathStr '/external/YALMIP']);
                addpath([pathStr '/external/YALMIP/demos']);
                addpath([pathStr '/external/YALMIP/extras']);
                addpath([pathStr '/external/YALMIP/modules']);
                addpath([pathStr '/external/YALMIP/modules/bilevel']);
                addpath([pathStr '/external/YALMIP/modules/global']);
                addpath([pathStr '/external/YALMIP/modules/moment']);
                addpath([pathStr '/external/YALMIP/modules/parametric']);
                addpath([pathStr '/external/YALMIP/modules/robust']);
                addpath([pathStr '/external/YALMIP/modules/sos']);
                addpath([pathStr '/external/YALMIP/operators']);
                addpath([pathStr '/external/YALMIP/solvers']);
                if verbose >= 1
                    disp('Adding embedded YALMIP to the path');
                end
            end
        elseif verbose >= 2
            disp('YALMIP is already in the path');
        end

        % Making sure SDPT3 is in the path and working
        SDPT3InPath = false;
        try
            [blk, Avec, C, b, X0, y0, Z0] = randsdp([2 2], [2 2], 2, 2);
            options.printlevel = 0;
            sdpt3(blk, Avec, C, b, options);
            SDPT3InPath = true;
        catch
        end
        if ~SDPT3InPath
            if exist([pathStr '/external/SDPT3/sdpt3.m']) ~= 2
                if verbose >= 2
                    warning(['The SDPT3 library was found neither in the path nor in the folder ', pathStr, '/external/SDPT3, some functionalities of the library might be disabled']);
                end
            else
                addpath([pathStr '/external/SDPT3']);
                if verbose >= 0
                    disp('Adding embedded SDPT3 to the path. Please run ''install_sdpt3'' to complete the setup of this solver');
                end
            end
        elseif verbose >= 2
            disp('SDPT3 is already in the path');
        end
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
