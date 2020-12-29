classdef Equivariant_forSubRep < replab.Equivariant
% Describes a subspace of an equivariant space induced by subspaces of source/target representations
%
% This equivariant space relates to a parent equivariant space as follows.
%
% - ``self.repR`` must be a `~+replab.SubRep` of ``parent.repR``,
% - ``self.repC`` must be a `~+replab.SubRep` of ``parent.repC``.

    properties (SetAccess = protected)
        parent % (`+replab.Equivariant`): Parent equivariant space
    end

    methods

        function self = Equivariant_forSubRep(parent, subRepR, subRepC, special)
        % Constructs an equivariant space similar to another equivariant space
        %
        % Args:
        %   parent (`+replab.Equivariant`): Equivariant space of ``subRepR.parent`` and ``subRepC.parent``
        %   subRepR (`+replab.SubRep`): A subrepresentation of ``parent.repR``
        %   subRepC (`+replab.SubRep`): A subrepresentation of ``parent.repC``
        %   special (charstring): Whether the equivariant subspace has special structure
            self@replab.Equivariant(subRepR, subRepC, special);
            assert(isa(parent, 'replab.Equivariant'));
            assert(isa(subRepR, 'replab.SubRep'));
            assert(isa(subRepC, 'replab.SubRep'));
            assert(subRepR.parent == parent.repR);
            assert(subRepC.parent == parent.repC);
            self.parent = parent;
        end

    end

    methods (Access = protected) % Implementations

        % Equivariant

        function X1 = project_exact(self, X)
            parentX = self.repR.injection('exact') * X * self.repC.projection('exact');
            X1 = self.repR.projection('exact') * self.parent.project(parentX) * self.repC.injection('exact');
        end

        function [X1, err] = project_double_sparse(self, X)
            parentX = self.repR.injection('double/sparse') * X * self.repC.projection('double/sparse');
            if nargout > 1
                [parentX1, eParent] = self.parent.project(parentX, 'double/sparse');
            else
                parentX1 = self.parent.project(parentX, 'double/sparse');
            end
            X1 = self.repR.projection('double/sparse') * parentX1 * self.repC.injection('double/sparse');
            if nargout > 1
                err = inf; % TODO
            end
        end

    end

end
