classdef MOcov < replab.init.ExternalDependency

    methods

        function self = MOcov
            self@replab.init.ExternalDependency('MOcov', 'MOcov/mocov.m');
        end

        function res = inPath(self)
            res = any(exist('mocov') == 2);
        end

        function res = works(self)
            res = ~isempty(mocov_get_absolute_path('.'));
        end

        function init(self, path)
            addpath(fullfile(path, 'MOcov'));
        end

    end

end
