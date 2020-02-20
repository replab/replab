classdef ForSubReps < replab.Equivariant
% Describes a subspace of an equivariant space induced by subspaces of source/target representations
%
% The child equivariant space (this object) relates to a parent equivariant space as follows.
%
% - ``child.repR`` must be a `~+replab.SubRep` of ``parent.repR``,
% - ``child.repC`` must be a `~+replab.SubRep` of ``parent.repC``.

    properties (SetAccess = protected)
        parent % (`+replab.Equivariant`): Parent equivariant space
    end

    methods

        function self = ForSubReps(repC, repR, special, parent)
            self@replab.Equivariant(repC, repR, special);
            assert(isa(repC, 'replab.SubRep'));
            assert(isa(repR, 'replab.SubRep'));
            assert(repC.parent == parent.repC);
            assert(repR.parent == parent.repR);
        end

        function [X1 err] = project(self, X)
            parentX = self.repR.H_internal * X * self.repC.F_internal;
            parentX1 = self.parent.project(X);
            X1 = full(self.repR.F_internal * parentX1 * self.repC.H_internal);
            err = NaN;
        end

        function [X err] = sampleWithError(self)
            parentX = self.parent.sampleWithError;
            X = full(self.repR.F_internal * parentX * self.repC.H_internal);
        end

    end

end
