classdef IsotypicEquivariant_trivial < replab.IsotypicEquivariant

    methods

        function self = IsotypicEquivariant_trivial(parent, repR, repC, special)
            self@replab.IsotypicEquivariant(parent, repR, repC, special, cell(1, 0), cell(1, 0));
        end

    end

    methods % Implementations

        % IsotypicEquivariant

        function b = isZero(self)
            b = true;
        end

        function X = reconstruct(self, M, type)
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            if strcmp(type, 'exact')
                X = replab.cyclotomic.zeros(self.repR.dimension, self.repC.dimension);
            else
                X = zeros(self.repR.dimension, self.repC.dimension);
            end
        end

        function [M, err] = projectAndFactor(self, X, type)
            M = cell(1, 0);
            err = 0;
        end

        function [M, err] = projectAndFactorFromParent(self, X, type)
            M = cell(1, 0);
            err = 0;
        end

    end

end
