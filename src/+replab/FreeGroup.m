classdef FreeGroup < replab.FPGroup
% Describes the free group

    methods (Static)

        function [F, varargout] = of(varargin)
            F = replab.FreeGroup(varargin, replab.FPGroup.newId);
            for i = 1:F.nGenerators
                varargout{i} = F.generator(i);
            end
        end

    end

    methods % superclass implementations

        function self = FreeGroup(names, id)
            self@replab.FPGroup(names, {}, id);
        end

        function x = sample(self)
            l = 10;
            x = randi([-self.n self.n], 1, l);
            x = self.reduce(x(x ~= 0));
        end

        function res = eqv(self, x, y)
            res = (isempty(x) && isempty(y)) || isequal(x, y);
        end

    end

end
