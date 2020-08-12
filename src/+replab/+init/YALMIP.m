classdef YALMIP < replab.init.ExternalDependency

    properties (Constant = true)
        subfolders = {'demos' 'extras' 'modules' 'modules/bilevel' 'modules/global' 'modules/moment' ...
                      'modules/parametric' 'modules/robust' 'modules/sos' 'operators' 'solvers'};
    end

    methods

        function self = YALMIP
            self@replab.init.ExternalDependency('YALMIP', 'yalmiptest.m');
        end

        function res = inPath(self)
            res = exist('yalmiptest') == 2;
        end

        function res = works(self)
            res = ~isempty(yalmip('version'));
        end

        function init(self, path)
            subfolders = self.subfolders;
            addpath(fullfile(path));
            for i = 1:length(subfolders)
                addpath(fullfile(path, subfolders{i}));
            end
        end
    end

end
