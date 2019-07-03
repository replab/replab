classdef Irreducible < replab.Str
% Describes the irreducible decomposition of a representation
    
    properties
        parent;     % Parent representation
        group;      % Group represented
        field;      % Field ('R' real or 'C' complex)
        components; % Isotypic components
    end

    methods

        function self = Irreducible(parent, components)
            self.parent = parent;
            self.group = parent.group;
            self.field = parent.field;
            self.components = components;
        end
        
        function n = nComponents(self)
            n = length(self.components);
        end
        
        function c = component(self, i)
            c = self.components{i};
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
