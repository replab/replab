classdef ComplexifiedRep < replab.Rep
% Representation derived by complexfying a real representation

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Real representation being transformed
    end

    methods

        function self = ComplexifiedRep(parent)
            assert(isequal(parent.field, 'R'), 'Only a real representation can be complexified.');
            args = cell(1, 0);
            if parent.inCache('isUnitary')
                args = horzcat(args, {'isUnitary' parent.isUnitary});
            end
            if parent.inCache('trivialDimension')
                args = horzcat(args, {'trivialDimension' parent.trivialDimension});
            end
            if parent.inCache('frobeniusSchurIndicator')
                args = horzcat(args, {'frobeniusSchurIndicator' parent.frobeniusSchurIndicator});
            end
            self@replab.Rep(parent.group, 'C', parent.dimension, args{:});
            self.parent = parent;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_complexifiedTrivial(self, options)
        % Complexified trivial is complex trivial
            if isa(self.parent, 'replab.rep.TrivialRep')
                res = replab.rep.TrivialRep(self.parent.group, 'C', self.parent.dimension);
            else
                res = [];
            end
        end

    end

    methods (Access = protected)

        % Rep

        function rep = computeDouble(self)
            rep = replab.rep.ComplexifiedRep(double(self.parent));
        end

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParents)
            r = replab.rep.ComplexifiedRep(newParents{1});
        end

        function rho = image_double_sparse(self, g)
            rho = self.parent.image(g, 'double/sparse');
        end

        function rho = image_exact(self, g)
            rho = self.parent.image(g, 'exact');
        end

        function e = computeErrorBound(self)
            e = self.parent.errorBound;
        end

        function c = computeConditionNumberEstimate(self)
            c = self.parent.conditionNumberEstimate;
        end

        function k = computeKernel(self)
            k = self.parent.kernel;
        end

        function rep = computeUnitarize(self)
            sr = self.parent.unitarize;
            rep = replab.SimilarRep(self, sr.A_internal, sr.Ainv_internal);
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = 'Complexification of representation';
        end

        % Rep

        function b = isExact(self)
            b = self.parent.isExact;
        end

    end

end
