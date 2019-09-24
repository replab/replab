classdef Irreducible < replab.Str
% Describes the irreducible decomposition of a representation
    
    properties
        parent % replab.Rep: Parent representation, must be unitary
        components % row cell array of replab.Isotypic: Isotypic components
    end

    methods

        function self = Irreducible(parent, components)
            assert(isequal(parent.isUnitary, true));
            self.parent = parent;
            self.components = components;
        end
        
        function r = asRep(self)
            r = replab.IrreducibleRep(self);
        end
        
        function U = U(self)
            U = zeros(0, self.parent.dimension);
            for i = 1:self.nComponents
                U = [U; self.component(i).rep.U];
            end
        end
        
        function r = asConjugateRep(self)
        % Returns the subrepresentation corresponding to this irreducible decomposition
            r = self.parent.leftConjugateUnitary(self.U);
        end

        function n = nComponents(self)
            n = length(self.components);
        end
        
        function c = component(self, i)
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
