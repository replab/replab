classdef vpi < replab.init.ExternalDependency
% Adds the VPI library to the path if it is not yet present

    methods

        function self = vpi
            self@replab.init.ExternalDependency('vpi', '@vpi/vpi.m');
            self.name = 'vpi';
        end

        function res = inPath(self)
            res = exist('vpi') == 2;
        end

        function res = works(self)
            res = false;
            try
                res = isequal(strtrim(num2str(vpi('1234567890098765432100123456789'))), '1234567890098765432100123456789');
            catch
            end
        end

        function init(self, folderName)
            addpath(folderName);
        end

    end

end
