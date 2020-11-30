classdef Equivariant_forSimilarRep < replab.Equivariant

    properties (SetAccess = protected)
        parent % (`+replab.Equivariant`): Parent equivariant space
    end

    properties (Access = protected)
        repC_internal % (`+replab.SimilarRep`): Representation of row space as a similar representation
        repR_internal % (`+replab.SimilarRep`): Representation of row space as a similar representation
    end

    methods

        function self = Equivariant_forSimilarRep(repC, repR, special)
            self@replab.Equivariant(repC, repR, special);
            repC_internal = repC.simplify;
            repR_internal = repR.simplify;
            if ~isa(repC_internal, 'replab.SimilarRep')
                repC_internal = replab.SimilarRep.identical(repC_internal);
            end
            if ~isa(repR_internal, 'replab.SimilarRep')
                repR_internal = replab.SimilarRep.identical(repR_internal);
            end
            switch special
              case 'hermitian'
                parent = repR_internal.parent.hermitianInvariant;
              case 'commutant'
                parent = repR_internal.parent.commutant;
              case 'trivialRows'
                parent = repC_internal.parent.trivialRowSpace;
              case 'trivialCols'
                parent = repR_internal.parent.trivialColSpace;
              otherwise
                parent = repC_internal.parent.equivariantTo(repR_internal.parent);
            end
            self.parent = parent;
            self.repR_internal = repR_internal;
            self.repC_internal = repC_internal;
        end

    end

    methods (Access = protected)

        function X = project_exact(self, X)
            parentX = self.repR_internal.Ainv('exact') * X * self.repC_internal.A('exact');
            X = self.repR_internal.A('exact') * self.parent.project(parentX) * self.repC_internal.Ainv('exact');
        end

        function [X1 eX] = project_double_sparse(self, X)
            parentX = self.repR_internal.Ainv('double/sparse') * X * self.repC_internal.A('double/sparse');
            if nargout > 1
                [parentX1 eParent] = self.parent.project(parentX, 'double/sparse');
            else
                parentX1 = self.parent.project(parentX, 'double/sparse');
            end
            X1 = self.repR_internal.A('double/sparse') * parentX1 * self.repC_internal.Ainv('double/sparse');
            if nargout > 1
                c = self.repR_internal.basisConditionNumberEstimate * self.repC_internal.basisConditionNumberEstimate;
                % error in the first conversion
                e1 = replab.numerical.norm2UpperBound(parentX) * self.repR_internal.A_Ainv_error * c;
                % error in the projection
                e2 = eParent * c;
                % error in the second conversion
                e3 = replab.numerical.norm2UpperBound(X1) * self.repC_internal.A_Ainv_error;
                eX = e1 + e2 + e3;
            end
        end

    end

end
