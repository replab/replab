classdef intlab < replab.init.Dependency
% Verifies that the INTLAB library is available and runs its setup if necessary

    methods

        function self = intlab
            self.name = 'intlab';
        end

        function res = inPath(self)
            res = any(exist('intval') == [2 6]);
        end

        function res = works(self)
            res = self.inPath;
        end

        function init(self, folderName)
            run(fullfile(path, 'startintlab.m'));
        end

    end

end
