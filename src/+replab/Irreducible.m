classdef Irreducible < replab.SubRep
% Describes the irreducible decomposition of a representation
%
% The mathematical background is based on Section 2.6 of Jean-Pierre Serre, "Linear representations of finite groups".
%
% An irreducible decomposition contains a sequence of isotypic components. When the irreducible decomposition is obtained
% through `+replab.Rep.decomposition` with an ``exact`` argument, the order of isotypic components matches the real or complex
% character table of the represented group, and thus the irreducible decomposition can contain empty isotypic components.
% Those empty components can be stripped using the `.squeeze` method.
%
% As the decomposition of a representation is an expensive operation in RepLAB, `.Irreducible` provides a pair of methods
% `.export` and `.import` that return a plain struct that can be saved into a ``.mat`` file under both MATLAB and Octave.

    properties (SetAccess = protected)
        components % (cell(1,\*) of `+replab.Isotypic`): Isotypic components
    end

    methods (Static)

        function I = import(parent, data)
        % Reconstruct an irreducible decomposition from exported data
        %
        % Args:
        %   parent (`.Rep`): Representation to decompose
        %   data (struct): Plain data resulting from `.export`
            assert(isequal(data.field, parent.field));
            assert(data.dimension == parent.dimension);
            assert(iscell(data.injection) == iscell(data.projection));
            if iscell(data.injection)
                injection = replab.cyclotomic(data.injection);
                projection = replab.cyclotomic(data.projection);
            else
                injection = data.injection;
                projection = data.projection;
            end
            n = length(data.components);
            shift = 0;
            components = cell(1, n);
            for i = 1:n
                compdata = data.components(i);
                comprange = shift+(1:compdata.irrepDimension*compdata.multiplicity);
                irreps = cell(1, compdata.multiplicity);
                for j = 1:compdata.multiplicity
                    subrange = shift+(1:compdata.irrepDimension);
                    shift = shift + compdata.irrepDimension;
                    subinj = injection(:, subrange);
                    subprj = projection(subrange, :);
                    irreps{j} = parent.subRep(subinj, 'projection', subprj, ...
                                              'frobeniusSchurIndicator', compdata.frobeniusSchurIndicator, ...
                                              'divisionAlgebraName', compdata.divisionAlgebraName, ...
                                              'isUnitary', compdata.isUnitary);
                end
                compprj = projection(comprange, :);
                components{i} = replab.Isotypic(parent, irreps, irreps{1}, compprj);
            end
            I = replab.Irreducible(parent, components);
        end

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

        function n = nComponents(self)
        % Returns the number of isotypic components in the decomposition
        %
        % Returns:
        %   integer: Number of isotypoic components
            n = length(self.components);
        end

        function data = export(self)
        % Returns a plain data structure with the information about this decomposition
        %
        % The returned data can be saved in a ``.mat`` file in both MATLAB and Octave and passed to the static
        % `.import` static method to reconstruct this irreducible decomposition.
        %
        % Returns:
        %   struct: Plain data structure
            nonEmptyComponents = cellfun(@(iso) iso.dimension > 0, self.components);
            assert(all(nonEmptyComponents), 'This irreducible decomposition must not contain any empty isotypic components. Call .squeeze.');
            if self.isExact
                inj = self.injection('exact').toCellOfStrings;
                proj = self.projection('exact').toCellOfStrings;
            else
                inj = self.injection;
                proj = self.projection;
            end
            n = self.nComponents;
            irrepDimension = cellfun(@(iso) iso.irrepDimension, self.components, 'uniform', 0);
            multiplicity = cellfun(@(iso) iso.multiplicity, self.components, 'uniform', 0);
            frobeniusSchurIndicator = cellfun(@(iso) iso.frobeniusSchurIndicator, self.components, 'uniform', 0);
            divisionAlgebraName = cellfun(@(iso) iso.divisionAlgebraName, self.components, 'uniform', 0);
            isUnitary = cellfun(@(iso) iso.isUnitary, self.components, 'uniform', 0);
            components = struct('irrepDimension', irrepDimension, 'multiplicity', multiplicity, ...
                                'frobeniusSchurIndicator', frobeniusSchurIndicator, 'divisionAlgebraName', divisionAlgebraName, ...
                                'isUnitary', isUnitary);
            data = struct('field', {self.field}, 'dimension', {self.dimension}, 'injection', {inj}, 'projection', {proj}, 'components', {components});
        end

        function [I, ind] = squeeze(self)
        % Returns an irreducible decomposition where empty isotypic components have been removed
        %
        % Useful to clean up the result of the exact decomposition, as the exact decomposition contains
        % isotypic components for all irreducible representations, even if they are not present in the
        % decomposed representation.
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> rep = S3.naturalRep;
        %   >>> dec = rep.decomposition('exact');
        %   >>> dec.nComponents
        %       3
        %   >>> dec.squeeze.nComponents
        %       2
        %
        % Returns
        % -------
        %   I: `.Irreducible`
        %     Irreducible decomposition with empty isotypic components removed
        %   ind: integer(1,\*)
        %     Indices of non-empty components
            mask = cellfun(@(iso) iso.dimension > 0, self.components);
            if all(mask)
                I = self;
            else
                I = replab.Irreducible(self.parent, self.components(mask));
            end
            ind = find(mask);
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

        function rep = conj(self)
            components = cellfun(@(iso) conj(iso), self.components, 'uniform', 0);
            rep = replab.Irreducible(conj(self.parent), components);
        end

        function rep = dual(self)
            components = cellfun(@(iso) dual(iso), self.components, 'uniform', 0);
            rep = replab.Irreducible(dual(self.parent), components);
        end

        % Rep: Equivariant spaces

        function a = antilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            a = self.cached(['antilinearInvariant_' type], @() self.irreducibleEquivariantFrom(conj(self),  'special', 'antilinear', 'type', type));
        end

        function b = bilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            b = self.cached(['bilinearInvariant_' type], @() self.irreducibleEquivariantTo(dual(self),  'special', 'bilinear', 'type', type));
        end

        function c = commutant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['commutant_' type], @() self.irreducibleEquivariantFrom(self,  'special', 'commutant', 'type', type));
        end

        function h = hermitianInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            h = self.cached(['hermitianInvariant_' type], @() self.irreducibleEquivariantTo(conj(dual(self)),  'special', 'hermitian', 'type', type));
        end

        function s = sesquilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            s = self.cached(['sesquilinearInvariant_' type], @() self.irreducibleEquivariantTo(conj(dual(self)),  'special', 'sesquilinear', 'type', type));
        end

        function b = symmetricInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            b = self.cached(['symmetricInvariant_' type], @() self.irreducibleEquivariantTo(dual(self),  'special', 'symmetric', 'type', type));
        end

        % SubRep

        function [res, anyBetter] = nice(self)
            components1 = cell(1, self.nComponents);
            anyBetter = false;
            for i = 1:self.nComponents
                [res, better] = self.components{i}.nice;
                components1{i} = res;
                anyBetter = anyBetter || better;
            end
            if anyBetter
                res = replab.Irreducible(self.parent, components1);
            else
                res = self;
            end
        end

        function irr = refine(self)
            components1 = cellfun(@(c) c.refine, self.components, 'uniform', 0);
            irr = replab.Irreducible(self.parent, components1);
        end

    end

end
