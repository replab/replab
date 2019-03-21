function replab_addpaths(includeEmbeddedSolver)
    % replab_addpaths([includeEmbeddedSolver])
    %
    % Adds directories that are needed to enable all functionalities of the
    % library to the search path.
    %
    % The optional parameter 'includeEmbeddedSolver' is 0 by default. Set
    % it to 1 to add the folders 'external/YALMIP/...' and 'externl/SeDuMi'
    % to the path (avoid this if you already have you installation of these
    % libraries)
    %
    % example: replab_addpaths
    
    if nargin < 1
        includeEmbeddedSolver = 0;
    end
    
    [pathStr, name, extension] = fileparts(which(mfilename));
    addpath(pathStr);
    addpath([pathStr '/external/vpi']);
    addpath([pathStr '/external/MOxUnit/MOxUnit']);
    addpath([pathStr '/external/MOxUnit/MOxUnit/util']);
    
    if includeEmbeddedSolver == 1
        % YALMIP
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
        
        % SDPT3
        addpath([pathStr '/external/SDPT3']);
    end
end
