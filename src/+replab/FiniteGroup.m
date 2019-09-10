classdef FiniteGroup < replab.CompactGroup
    
    properties (SetAccess = protected)
        generators % row cell array of group elements: Group generators
        order % vpi: Order of the group, equal to the number of distinct group elements
        niceMonomorphism % function_handle: Injective group homomorphism from this group into a permutation group
    end
    
    properties (Access = protected)
        randomBag_ % Bag of elements used for random sampling
    end
    
    methods
        
        function R = randomBag(self)
            if isequal(self.randomBag_, [])
                self.randomBag_ = replab.RandomBag(self, self.generators);
            end
            R = self.randomBag_;
        end
        
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Group(self);
            names{1, end+1} = 'generators';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nGenerators
                names{1, end+1} = sprintf('generator(%d)', i);
                values{1, end+1} = self.generator(i);
            end
        end
        
        function n = nGenerators(self)
        % Returns the number of group generators
        %
        % Returns:
        %   integer: Number of group generators
            n = length(self.generators);
        end
        
        function p = generator(self, i) 
        % Returns the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: i-th group generator
            p = self.generators{i};
        end
        
        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: Inverse of the i-th group generator
            p = self.inverse(self.generators{i});
        end

        function b = isTrivial(self)
        % Tests whether this group is trivial
        %
        % Returns:
        %   logical: True if this group is trivial (i.e. has only one element)
            b = self.nGenerators == 0;
        end
        
        function o = order(self)
        % Returns the group order, computing it if necessary
        %
        % Returns:
        %   vpi: Order of this group
            error('Not implemented');
        end
        
        function e = elements(self)
        % Returns an enumeration of the group elements
        %
        % Returns:
        %   replab.Enumerator: A space-efficient enumeration of the group elements
            error('Not implemented');
        end
        
        function D = decomposition(self)
        % Returns a decomposition of this group as a product of sets
        %
        % Returns:
        %   replab.FiniteGroupDecomposition: The group decomposition
            assert(self.order < 1e6, 'Default decomposition is available only for small groups');
            O = double(self.order);
            C = self.elements.toCell;
            idIndex = self.elements.find(self.identity);
            C = C([idIndex setdiff(1:O, idIndex)]);
            D = replab.FiniteGroupDecomposition.trivial(self, C);
        end

        function g = sampleUniformly(self)
            g = self.elements.sample;
        end

        function rep = leftRegularRep(self)
            o = self.order;
            assert(o < 1e6);
            o = double(o);
            perms = cell(1, self.nGenerators);
            E = self.elements;
            for i = 1:self.nGenerators
                g = self.generator(i);
                img = zeros(1, o);
                for j = 1:o
                    img(j) = double(E.find(self.compose(g, E.at(j))));
                end
                perms{i} = img;
            end
            rep = self.permutationRep(o, perms);
        end

    end

end
