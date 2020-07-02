function res = initYALMIP(verbose)
% Adds the YALMIP library to the path if it is not yet present
%
% Args:
%   verbose ({0, 1, 2}): Controls the display level
%
% Returns:
%   logical: True if the library is available

    basePath = replab.globals.replabPath;
    res = false;
    try
        yalmip('version');
        res = true;
    catch
    end
    if ~res
        if exist([basePath '/external/YALMIP/yalmiptest.m']) ~= 2
            warning(['The YALMIP library was found neither in the path nor in the folder ', basePath, '/external/YALMIP.', char(10), ...
                'Some functionalities of the library might be disabled.']);
        else
            addpath([basePath '/external/YALMIP']);
            addpath([basePath '/external/YALMIP/demos']);
            addpath([basePath '/external/YALMIP/extras']);
            addpath([basePath '/external/YALMIP/modules']);
            addpath([basePath '/external/YALMIP/modules/bilevel']);
            addpath([basePath '/external/YALMIP/modules/global']);
            addpath([basePath '/external/YALMIP/modules/moment']);
            addpath([basePath '/external/YALMIP/modules/parametric']);
            addpath([basePath '/external/YALMIP/modules/robust']);
            addpath([basePath '/external/YALMIP/modules/sos']);
            addpath([basePath '/external/YALMIP/operators']);
            addpath([basePath '/external/YALMIP/solvers']);
            if verbose >= 1
                disp('Adding embedded YALMIP to the path');
            end
            res = true;
        end
    elseif verbose >= 2
        disp('YALMIP is already in the path');
    end
end
