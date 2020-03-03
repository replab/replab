classdef Irreducible < replab.SubRep
% Describes the irreducible decomposition of a representation
%
% For the background, see Section 2.6 of Jean-Pierre Serre, Linear representations of finite groups
%
% The irreducible decomposition of ``parent`` contains isotypic components in the cell vector ``components``.
% Each isotypic component corresponds to a set of equivalent irreducible representations.

    properties
        components % (cell(1,\*) of `+replab.Isotypic`): Isotypic components
    end

    methods

        function self = Irreducible(parent, components)
            n = length(components);
            Bs = cell(1, n);
            Es = cell(n, 1);
            for i = 1:n
                c = components{i};
                assert(isa(c, 'replab.Isotypic'));
                Bs{1,i} = c.B_internal;
                Es{i,1} = c.E_internal;
            end
            B_internal = [Bs{:}];
            E_internal = vertcat(Es{:});
            self = self@replab.SubRep(parent, B_internal, E_internal);
            self.components = components;
        end

        function r = asSimilarRep(self)
        % Returns the block-diagonal representation corresponding to the decomposition
        %
        % Up to the change of basis matrix `U`, it corresponds to the representation ``parent``.
        % Indeed, we have ``self.image(g) = U * self.parent.image(g) * U'``.
        %
        % The returned representation is the parent representation with an explicit change of basis, so it does
        % not look as clean as ``self``. For efficiency and numerical stability, use ``self``.
        %
        % Returns:
        %   `+replab.Rep`: The block-diagonal representation as a representation similar to this rep. parent
            r = self.parent.similarRep(self.B_internal, self.E_internal);
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

        %% Str methods

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

        %% Rep methods

        function rho = image_internal(self, g)
            blocks = cellfun(@(iso) iso.image_internal(g), self.components, 'uniform', 0);
            % Construct the blocks in the block diagonal image
            rho = blkdiag(blocks{:});
        end

        function c = commutant(self)
            if isempty(self.commutant_)
                self.commutant_ = replab.IrreducibleCommutant(self);
            end
            c = self.commutant_;
        end

        %% SubRep methods

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
