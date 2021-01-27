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

    methods % Simplfication rules

        % Rep

        function res = rewriteTerm_distributeOverDirectSum(self, options)
        % conj(.), and dual(.) distribute over a direct sum
            if isa(self.parent, 'replab.rep.DirectSumRep')
                newFactors = cellfun(@(f) replab.rep.DerivedRep(f, self.conjugate, self.inverse, self.transpose), ...
                                  self.parent.factors, 'uniform', 0);
                res = replab.rep.DirectSumRep(self.parent.group, self.parent.field, newFactors);
            else
                res = [];
            end
        end

        function res = rewriteTerm_distributeOverTensorProduct(self, options)
        % conj(.), and dual(.) distribute over a tensor product
            if isa(self.parent, 'replab.rep.TensorRep')
                newFactors = cellfun(@(f) replab.rep.DerivedRep(f, self.conjugate, self.inverse, self.transpose), ...
                                     self.parent.factors, 'uniform', 0);
                res = replab.rep.TensorRep(self.parent.group, self.parent.field, newFactors);
            else
                res = [];
            end
        end

        function res = rewriteTerm_noEffect(self, options)
        % Remove this DerivedRep if it has no effect
            if ~self.conjugate && ~self.inverse && ~self.transpose
                res = self.parent;
            else
                res = [];
            end
        end

        function res = rewriteTerm_conjugateReal(self, options)
        % Remove the conjugation part if it applies to a real representation
            if self.overR && self.conjugate
                res = replab.rep.DerivedRep(self.parent, false, self.inverse, self.transpose);
            else
                res = [];
            end
        end

        function res = rewriteTerm_dualIsConjugateWhenUnitary(self, options)
        % rho_{g^-1}.' = conj(rho) if rho is unitary
            if self.parent.knownUnitary && self.inverse && self.transpose
                res = replab.rep.DerivedRep(self.parent, ~self.conjugate, false, false);
            else
                res = [];
            end
        end

        function res = rewriteTerm_parentIsTrivialRep(self, options)
        % We have no effect when applied on the trivial representation
            if isa(self.parent, 'replab.rep.TrivialRep')
                res = self.parent;
            else
                res = [];
            end
        end

        function res = rewriteTerm_parentIsComplexifiedRep(self, options)
        % Move DerivedRep inside ComplexifiedRep where it can be simplified
            if isa(self.parent, 'replab.rep.ComplexifiedRep')
                res = replab.rep.ComplexifiedRep(replab.rep.DerivedRep(self.parent.parent, ...
                                                                  self.conjugate, self.inverse, self.transpose));
            else
                res = [];
            end
        end

        function res = rewriteTerm_parentIsDerivedRep(self, options)
        % Collapse two successive DerivedRep
            if isa(self.parent, 'replab.rep.DerivedRep')
                res = replab.rep.DerivedRep(self.parent.parent, ...
                                            xor(self.parent.conjugate, self.conjugate), ...
                                            xor(self.parent.inverse, self.inverse), ...
                                            xor(self.parent.transpose, self.transpose));
            else
                res = [];
            end
        end

        function res = rewriteTerm_parentIsSimilarRep(self, options)
            if isa(self.parent, 'replab.SimilarRep')
                A = self.parent.A_internal;
                Ainv = self.parent.Ainv_internal;
                if self.conjugate
                    A = conj(A);
                    Ainv = conj(Ainv);
                end
                swap = false;
                if self.transpose
                    A0 = A;
                    Ainv0 = Ainv;
                    A = Ainv0.';
                    Ainv = A0.';
                end
                res = replab.SimilarRep(replab.rep.DerivedRep(self.parent.parent, ...
                                                              self.conjugate, self.inverse, self.transpose), ...
                                        A, Ainv);
            else
                res = [];
            end
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rep = computeDouble(self)
            rep = replab.rep.DerivedRep(double(self.parent), self.conjugate, self.inverse, self.transpose);
        end

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParents)
            r = replab.rep.DerivedRep(newParents{1}, self.conjugate, self.inverse, self.transpose);
        end

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

        function b = computeIsUnitary(self)
            b = self.parent.isUnitary;
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

        function p = invariantBlocks(self)
            p = self.parent.invariantBlocks; % the block structure does not change
        end

        function res = dual(self)
            res = replab.rep.DerivedRep(self.parent, self.conjugate, ~self.inverse, ~self.transpose);
        end

        function res = conj(self)
            res = replab.rep.DerivedRep(self.parent, ~self.conjugate, self.inverse, self.transpose);
        end

    end

end
