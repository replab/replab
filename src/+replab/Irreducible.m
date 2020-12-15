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
            injections = cell(1, n);
            projections = cell(1, n);
            for i = 1:n
                c = components{i};
                assert(isa(c, 'replab.Isotypic'));
                injections{i} = c.injection_internal;
                projections{i} = c.projection_internal;
            end
            injection = horzcat(injections{:});
            projection = vertcat(projections{:});
            self = self@replab.SubRep(parent, injection, projection);
            self.components = components;
            if length(components) == 1 && components{1}.multiplicity == 1
                self.cache('isIrreducible', true, '==');
            end
            %if isequal(self.basis, eye(self.dimension))
            %    replab.rep.copyProperties(self.components{1}.irreps{1}, self);
            %end
        end

        function r = asSimilarRep(self)
        % Returns the block-diagonal similar representation corresponding to the decomposition
        %
        % The returned representation is the parent representation with an explicit change of basis, so it does
        % not look as clean as ``self``. For efficiency and numerical stability, use ``self``.
        %
        % Returns:
        %   `+replab.Rep`: The block-diagonal representation as a representation similar to this rep. parent
            r = self.parent.similarRep(self.projection_internal, self.injection_internal);
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

% $$$         % Rep
% $$$
% $$$         function rho = image_internal(self, g)
% $$$             blocks = cellfun(@(iso) iso.image_internal(g), self.components, 'uniform', 0);
% $$$             % Construct the blocks in the block diagonal image
% $$$             rho = blkdiag(blocks{:});
% $$$         end

% $$$         function c = computeCommutant(self)
% $$$             c = replab.IrreducibleCommutant(self);
% $$$         end

% $$$         % SubRep
% $$$
% $$$         function irr = nice(self)
% $$$             components1 = cellfun(@(c) c.nice, self.components, 'uniform', 0);
% $$$             irr = replab.Irreducible(self.parent, components1);
% $$$         end
% $$$
% $$$         function irr = refine(self)
% $$$             components1 = cellfun(@(c) c.refine, self.components, 'uniform', 0);
% $$$             irr = replab.Irreducible(self.parent, components1);
% $$$         end

    end

end
