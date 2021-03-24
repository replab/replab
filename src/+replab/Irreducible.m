classdef Irreducible < replab.SubRep
% Describes the irreducible decomposition of a representation
%
% For the background, see Section 2.6 of Jean-Pierre Serre, Linear representations of finite groups
%
% The irreducible decomposition of ``parent`` contains isotypic components in the cell vector ``components``.
% Each isotypic component corresponds to a set of equivalent irreducible representations.

    properties (SetAccess = protected)
        components % (cell(1,\*) of `+replab.Isotypic`): Isotypic components
    end

    methods

        function self = Irreducible(parent, components)
            n = length(components);
            injections = cellfun(@(c) c.injection_internal, components, 'uniform', 0);
            projections = cellfun(@(c) c.projection_internal, components, 'uniform', 0);
            isUnitary = cellfun(@(c) c.isUnitary, components);
            assert(all(cellfun(@(c) isa(c, 'replab.Isotypic'), components)), 'All components must be isotypic');
            injection = horzcat(injections{:});
            projection = vertcat(projections{:});
            isIrreducible = length(components) == 1 && components{1}.multiplicity == 1;
            self = self@replab.SubRep(parent, injection, projection, 'isUnitary', all(isUnitary), 'isIrreducible', isIrreducible);
            self.components = components;
        end

        function r = toSubRep(self)
        % Returns the block-diagonal similar representation corresponding to this decomposition
        %
        % Both this `.Irreducible` object and the returned representation will be instances of `.SubRep`.
        % This method helps us test the implementation of `.Irreducible` as the returned `.SubRep` loses
        % its particular structure.
        %
        % The returned representation is the parent representation with an explicit change of basis, so it does
        % not look as clean as this `.Irreducible`.
        %
        % Returns:
        %   `+replab.SubRep`: The block-diagonal representation as a representation similar to `.parent`
            r = self.parent.similarRep(self.projection_internal, 'inverse', self.injection_internal);
        end

        function n = nComponents(self)
        % Returns the number of isotypic components in the decomposition
        %
        % Returns:
        %   integer: Number of isotypoic components
            n = length(self.components);
        end

        function c = component(self, i)
        % Returns a particular isotypic component in the decomposition
        %
        % Args:
        %   i (logical): Index of the isotypic component
        %
        % Returns:
        %   `+replab.Isotypic`: The ``i``-th isotypic component
            c = self.components{i};
        end

        function r = irrep(self, i, j)
        % Returns a subrepresentation in the irreducible decomposition
        %
        % Args:
        %   i (integer): Index of the isotypic component
        %   j (integer, optional): Index of the copy in the ``i``-th isotypic component
        %                          Default value is ``1``.
        %
        % Returns:
        %   `+replab.SubRep`: An irreducible subrepresentation
            if nargin < 3
                j = 1;
            end
            r = self.component(i).irrep(j);
        end

    end

    methods % Equivariant spaces

        function E = irreducibleEquivariantFrom(self, repC, varargin)
        % Returns the space of equivariant linear maps from another irreducible decomposition to this irreducible decomposition
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * repC.image(g) = self.image(g) * X``
        %
        % Both irreducible decompositions must be harmonized.
        %
        % Args:
        %   repC (`+replab.Irreducible`): Irreducible decomposition, representation on the source/column space
        %
        % Keyword Args:
        %   special (charstring, optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.IrreducibleEquivariant` or ``[]``: The equivariant vector space, or ``[]`` if the space has dimension zero or contains only the zero matrix
            E = replab.IrreducibleEquivariant.make(self, repC, varargin{:});
        end

        function E = irreducibleEquivariantTo(self, repR, varargin)
        % Returns the space of equivariant linear maps from this irreducible decomposition to another irreducible decomposition
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * self.image(g) = repR.image(g) * X``
        %
        % Both irreducible decompositions must be harmonized.
        %
        % Args:
        %   repR (`+replab.Irreducible`): Irreducible decomposition, representation on the target/row space
        %
        % Keyword Args:
        %   special ('commutant', 'hermitian', 'trivialRows', 'trivialCols' or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.IrreducibleEquivariant` or ``[]``: The equivariant vector space, or ``[]`` if the space has dimension zero or contains only the zero matrix
            E = replab.IrreducibleEquivariant.make(repR, self, varargin{:});
        end

    end

    methods (Access = protected) % Implementations


        function rho = image_double_sparse(self, g)
            blocks = cellfun(@(c) c.image(g, 'double/sparse'), self.components, 'uniform', 0);
            rho = blkdiag(blocks{:});
        end

        function rho = image_exact(self, g)
            blocks = cellfun(@(c) c.image(g, 'exact'), self.components, 'uniform', 0);
            rho = blkdiag(blocks{:});
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.SubRep(self);
            names{1, end+1} = 'components';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.SubRep(self);
            for i = 1:self.nComponents
                names{1, end+1} = sprintf('component(%d)', i);
                values{1, end+1} = self.component(i);
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.IrreducibleLaws(self);
        end

        % Rep

        % TODO: add other equivariant spaces
        function c = commutant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['commutant_' type], @() self.irreducibleEquivariantFrom(self,  'special', 'commutant', 'type', type));
        end

        function c = sesquilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['sesquilinearInvariant_' type], @() self.irreducibleEquivariantFrom(self,  'special', 'sesquilinear', 'type', type));
        end

         % SubRep

        function irr = nice(self)
            components1 = cellfun(@(c) c.nice, self.components, 'uniform', 0);
            irr = replab.Irreducible(self.parent, components1);
        end

        function irr = refine(self)
            components1 = cellfun(@(c) c.refine, self.components, 'uniform', 0);
            irr = replab.Irreducible(self.parent, components1);
        end

    end

end
