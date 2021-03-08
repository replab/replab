classdef SemidirectProductRep < replab.Rep
% Representation of a semidirect product using a representation of the factors
%
% Both factor representations act on the same vector space, and they are applied sequentially.

    properties (SetAccess = protected)
        Hrep % (`+replab.Rep`): Representation of the acting group
        Nrep % (`+replab.Rep`): Representation of the group acted upon
    end

    methods

        function self = SemidirectProductRep(group, Hrep, Nrep)
        % Constructs a representation of a semidirect product from a representation of the factors
        %
        % Args:
        %   group (`+replab.SemidirectProductGroup): Semidirect product group
        %   Hrep (`+replab.Rep`): Representation of ``group.H``
        %   Nrep (`+replab.Rep`): Representation of ``group.N``
            assert(isa(group, 'replab.SemidirectProductGroup'));
            assert(isa(Hrep, 'replab.Rep') && Hrep.group == group.H);
            assert(isa(Nrep, 'replab.Rep') && Nrep.group == group.N);
            assert(Hrep.field == Nrep.field);
            assert(Hrep.dimension == Nrep.dimension);
            isUnitary = Hrep.isUnitary && Nrep.isUnitary;
            self@replab.Rep(group, Hrep.field, Hrep.dimension, 'isUnitary', isUnitary);
            self.Hrep = Hrep;
            self.Nrep = Nrep;
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function c = decomposeTerm(self)
            c = {self.Hrep self.Nrep};
        end

        function r = composeTerm(self, newFactors)
            r = replab.rep.SemidirectProductRep(self.group, newFactors{1}, newFactors{2});
        end

        function rho = image_exact(self, g)
            rho = self.Hrep.image(g{1}, 'exact') * self.Nrep.image(g{2}, 'exact');
        end

        function rho = image_double_sparse(self, g)
            rho = self.Hrep.image(g{1}, 'double/sparse') * self.Nrep.image(g{2}, 'double/sparse');
        end

        function e = computeErrorBound(self)
            e1 = self.Hrep.errorBound;
            e2 = self.Nrep.errorBound;
            c1 = self.Hrep.conditionNumberEstimate;
            c2 = self.Nrep.conditionNumberEstimate;
            e = e1*c2 + c1*e2;
        end

        function c = computeConditionNumberEstimate(self)
            c = self.Hrep.conditionNumberEstimate * self.Nrep.conditionNumberEstimate;
        end

        function M = matrixRowAction_double_sparse(self, g, M)
            M = self.Hrep.matrixRowAction(g{1}, self.Nrep.matrixRowAction(g{2}, M));
        end

        function M = matrixColAction_double_sparse(self, g, M)
            M = self.Hrep.matrixColAction(g{1}, self.Nrep.matrixColAction(g{2}, M));
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = self.Hrep.isExact && self.Nrep.isExact;
        end

        function p = invariantBlocks(self)
            p = self.Hrep.invariantBlocks.join(self.Nrep.invariantBlocks);
        end

        function b = hasTorusImage(self)
            if ~self.Hrep.hasTorusImage || ~self.Nrep.hasTorusImage
                b = false;
                return
            end
            b = self.group.H.maximalTorusDimension == 0;
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            assert(self.hasTorusImage);
            [torusMap, torusInjection, torusProjection] = self.Nrep.torusImage;
        end

    end

end
