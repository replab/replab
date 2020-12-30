classdef Equivariant_forSimilarRep < replab.Equivariant
% Describes an equivariant space induced by similarity transformations on the source and target space.
%
% This equivariant space (this object) relates to a parent equivariant space as follows.
%
% - ``self.repR`` must be a `~+replab.SimilarRep` of ``parent.repR``,
% - ``self.repC`` must be a `~+replab.SimilarRep` of ``parent.repC``.

    properties (SetAccess = protected)
        parent % (`+replab.Equivariant`): Parent equivariant space
    end

    methods

        function self = Equivariant_forSimilarRep(parent, similarRepR, similarRepC, special)
        % Constructs an equivariant space similar to another equivariant space
        %
        % Args:
        %   parent (`+replab.Equivariant`): Equivariant space of ``similarRepR.parent`` and ``similarRepC.parent``
        %   similarRepR (`+replab.SimilarRep`): A similar representation to ``parent.repR``
        %   similarRepC (`+replab.SimilarRep`): A similar representation to ``parent.repC``
        %   special (charstring): Whether the equivariant subspace has special structure
             assert(isa(parent, 'replab.Equivariant'));
             assert(isa(similarRepR, 'replab.SimilarRep'));
             assert(isa(similarRepC, 'replab.SimilarRep'));
             %assert(similarRepR.parent == parent.repR);
             %assert(similarRepC.parent == parent.repC);
             self@replab.Equivariant(similarRepR, similarRepC, special);
             self.parent = parent;
         end

    end

    methods (Access = protected)

        function X1 = project_exact(self, X)
            parentX = self.repR.Ainv('exact') * X * self.repC.A('exact');
            X1 = self.repR.A('exact') * self.parent.project(parentX) * self.repC.Ainv('exact');
        end

        function [X1, err] = project_double_sparse(self, X)
            parentX = self.repR.Ainv('double/sparse') * X * self.repC.A('double/sparse');
            if nargout > 1
                [parentX1, eParent] = self.parent.project(parentX, 'double/sparse');
            else
                parentX1 = self.parent.project(parentX, 'double/sparse');
            end
            X1 = self.repR.A('double/sparse') * parentX1 * self.repC.Ainv('double/sparse');
            if nargout > 1
                c = self.repR.basisConditionNumberEstimate * self.repC.basisConditionNumberEstimate;
                % error in the first conversion
                e1 = replab.numerical.norm2UpperBound(parentX) * self.repR.A_Ainv_error * c;
                % error in the projection
                e2 = eParent * c;
                % error in the second conversion
                e3 = replab.numerical.norm2UpperBound(X1) * self.repC.A_Ainv_error;
                err = e1 + e2 + e3;
            end
        end

    end

end
