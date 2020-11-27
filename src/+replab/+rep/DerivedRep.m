classdef DerivedRep < replab.Rep
% Representation derived by the means of complex conjugate, inverse and/or transpose

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation being transformed
        conjugate % (logical): Whether to take the complex conjugate
        inverse % (logical): Whether to take the inverse of the group element
        transpose % (logical): Whether to transpose the map
    end

    methods

        function self = DerivedRep(parent, conjugate, inverse, transpose)
            assert(inverse == transpose, ...
                   'Cannot use inverse & transpose independently, as our representations are left modules');
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
            self@replab.Rep(parent.group, parent.field, parent.dimension, args{:});
            self.parent = parent;
            self.conjugate = conjugate;
            self.inverse = inverse;
            self.transpose = transpose;
        end

    end

    methods (Access = protected)

        % Rep

        function rho = image_double_sparse(self, g)
            if self.inverse
                g = self.group.inverse(g);
            end
            rho = self.parent.image(g, 'double/sparse');
            if self.conjugate
                rho = conj(rho);
            end
            if self.transpose
                rho = rho.';
            end
        end

        function rho = image_exact(self, g)
            if self.inverse
                g = self.group.inverse(g);
            end
            rho = self.parent.image(g, 'exact');
            if self.conjugate
                rho = conj(rho);
            end
            if self.transpose
                rho = rho.';
            end
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

        %TODO: optimize this
        %function rep = computeUnitarize(self)
        %    sr = self.parent.unitarize;
        %    rep = replab.SimilarRep(self, sr.A_internal, sr.Ainv_internal);
        %end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        % Rep

        function b = isExact(self)
            b = self.parent.isExact;
        end

    end

end
