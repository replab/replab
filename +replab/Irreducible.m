classdef Irreducible < replab.Str
% Describes the irreducible decomposition of a representation
    
    properties
        parent;     % Parent representation
        components; % Isotypic components
    end

    methods

        function self = Irreducible(parent, components)
            self.parent = parent;
            self.components = components;
        end
        
        function r = rep(self)
        % Returns the subrepresentation corresponding to this isotypic component
            U = zeros(self.parent.dimension, 0);
            for i = 1:self.nComponents
                U = [U self.component(i).rep.U];
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
        
        function I = recoverRational(self)
            components1 = cellfun(@(x) x.recoverRational, self.components, 'uniform', 0);
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
