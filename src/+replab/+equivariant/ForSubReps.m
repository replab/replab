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
            self.parent = parent;
        end

        function [X1 err] = project(self, X)
            parentX = self.repR.B_internal * X * self.repC.E_internal;
            parentX1 = self.parent.project(parentX);
            X1 = full(self.repR.E_internal * parentX1 * self.repC.B_internal);
            err = NaN;
        end

        function [X err] = sampleWithError(self)
            [parentX parentErr] = self.parent.sampleWithError;
            X = full(self.repR.E_internal * parentX * self.repC.B_internal);
            err = NaN; % TODO: include error of the basis
        end

        function [X err] = sampleInContext(self, context, ind)
            [parentX parentErr] = self.parent.sampleInContext(context, ind);
            X = full(self.repR.E_internal * parentX * self.repC.B_internal);
            err = NaN; % TODO: include error of the basis
        end

    end

    methods (Static)

        function e = make(repC, repR, special)
            if ~isa(repR, 'replab.SubRep') || ~ismember(special, {'commutant', 'hermitian', 'trivial'})
                e = replab.DispatchNext('Can only handle commutant or hermitian invariant spaces of subrepresentations.');
                return
            end
            switch special
              case 'commutant'
                parentRep = repC.parent;
                e = parentRep.commutant.subEquivariant(repC, repR, special);
              case 'hermitian'
                % parentR is the representation on which .hermitianInverse has been called
                parentH = repR.parent.hermitianInvariant;
                parentR = parentH.repR;
                parentC = parentH.repC;
                assert(parentR == repR.parent);
                repC1 = replab.SubRep(parentC, repR.E_internal', repR.B_internal');
                e = parentH.subEquivariant(repC1, repR, special);
              case 'trivial'
                % repC is trivial, repR.parent is non trivial
                assert(isa(repC, 'replab.rep.TrivialRep'));
                d = repC.dimension;
                assert(d == repR.dimension);
                parentT = repR.parent.trivialSpace;
                dParent = parentT.repC.dimension;
                H = sparse(1:d, 1:d, ones(1, d), dParent, d);
                F = H';
                repC1 = replab.SubRep(parentT.repC, H, F);
                e = parentT.subEquivariant(repC1, repR, special);
            end
        end

    end

end
