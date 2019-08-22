classdef Irreducible < replab.Str
% Describes the irreducible decomposition of a representation
%
% It corresponds to the canonical decomposition of Section 2.6 in :cite:`Serre1977`.
    properties
        parent % Representation being decomposed
        components % (row cell vector of :class:`+replab.Isotypic`): Isotypic components
    end

    methods

        function self = Irreducible(parent, components)
        % Constructor
        %
        % Parameters:
        %   parent (:class:`+replab.Rep`): Representation being decomposed
        %   components (row cell vector of :class:`+replab.Isotypic`): Isotypic components
            
            self.parent = parent;
            self.components = components;
        end
        
        function r = rep(self)
        % Returns the decomposed representation in the basis that expresses the decomposition
        %
        % Returns: 
        % Returns the subrepresentation corresponding to this isotypic component
            U = zeros(0, self.parent.dimension);
            for i = 1:self.nComponents
                U = [U; self.component(i).rep.U];
            end
            r = self.parent.leftConjugate(U);
            % TODO: preserve rational bases
        end

        function n = nComponents(self)
            n = length(self.components);
        end
        
        function c = component(self, i)
            c = self.components{i};
        end
        
        function I = nice(self)
            components1 = cellfun(@(x) x.nice, self.components, 'uniform', 0);
            I = replab.Irreducible(self.parent, components1);
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'components';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            for i = 1:self.nComponents
                names{1, end+1} = sprintf('component(%d)', i);
                values{1, end+1} = self.component(i);
            end
        end
        
    end

end
