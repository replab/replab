function res = initVPI
% Adds the VPI library to the path if it is not yet present
%
% Returns:
%   logical: True if the library is available

    verbose = replab.globals.verboseInit;
    basePath = replab.globals.replabPath;
    res = false;
    try
        res = isequal(strtrim(num2str(vpi('1234567890098765432100123456789'))), '1234567890098765432100123456789');
    catch
    end
    if ~res
        if exist([basePath '/external/vpi/@vpi/vpi.m']) ~= 2
            error(['The VPI library was not found in the folder ', basePath, '/external/vpi']);
        else
            addpath([basePath '/external/vpi']);
            if verbose >= 1
                disp('Adding VPI to the path');
            end
            res = true;
        end
    elseif verbose >= 2
        disp('VPI is already in the path');
    end
end
