classdef Irreducible < replab.SubRep
% Describes the irreducible decomposition of a representation
%
% For the background, see Section 2.6 of Jean-Pierre Serre, Linear representations of finite groups
%
% The irreducible decomposition of `parent` contains isotypic components in the cell vector `components`.
% Each isotypic component corresponds to a set of equivalent irreducible representations expressed in the same basis.
    
    properties
        components % row cell array of replab.Isotypic: Isotypic components
    end
    
    methods

        function self = Irreducible(parent, components)
            assert(isequal(parent.isUnitary, true));
            for i = 1:length(components)
                c = components{i};
                assert(isa(c, 'replab.Isotypic'));
            end
            Us = cellfun(@(iso) iso.U, components, 'uniform', 0);
            U = vertcat(Us{:});
            nbs = cellfun(@(iso) iso.niceBasis, components, 'uniform', 0);
            niceBasis = replab.NiceBasis.vertcat(nbs);
            self = self@replab.SubRep(parent, U, niceBasis);
            self.components = components;
        end
        
        function r = asConjugateRep(self)
        % Returns the block-diagonal representation corresponding to the decomposition
        %
        % Up to the change of basis matrix `self.U`, it corresponds to the representation `parent`.
        % Indeed, we have ``self.asRep.image(g) = U * self.parent.image(g) * U'``.
        %
        % The returned representation is a conjugate of the parent representation, so it does
        % not look as clean as `self.asRep`. For efficiency and numerical stability, use `self.asRep`.
        %
        % Returns:
        %   replab.Rep: The block-diagonal representation as a left conjugate representation
            r = self.parent.leftConjugateUnitary(self.U);
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
        %   replab.Isotypic: The `i`-th isotypic component
            c = self.components{i};
        end
        
        function r = irrep(self, i, j)
        % Returns a subrepresentation in the irreducible decomposition
        %
        % Args:
        %   i (integer): Index of the isotypic component
        %   j (integer, optional): Index of the copy in the `i`-th isotypic component
        %                          Default value is `1`.
        %
        % Returns:
        %   replab.SubRep: An irreducible subrepresentation
            if nargin < 3
                j = 1;
            end
            r = self.component(i).copy(j);
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
        
        function rho = image(self, g)
            blocks = cellfun(@(iso) iso.image(g), self.components, 'uniform', 0);
            % Construct the blocks in the block diagonal image
            rho = blkdiag(blocks{:});
        end

        function c = commutant(self)
            if isempty(self.commutant_)
                self.commutant_ = replab.IrreducibleCommutant(self);
            end
            c = self.commutant_;
        end
        
    end

end
