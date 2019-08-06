classdef FiniteGroup < replab.FinitelyGeneratedGroup
    
    methods
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.FinitelyGeneratedGroup(self);
            if self.knownOrder
                names{1, end+1} = 'order';
                values{1, end+1} = self.order;
            end
        end
        
        function b = contains(self, g)
        % Returns true when this group contains the element "g"
        % and false otherwise
            f = self.containsFun;
            b = f(g);
        end
        
        function b = knownOrder(self)
        % Returns true if this group order has been computed
            f = self.knownOrderFun;
            b = f();
        end
        
        function o = order(self)
        % Returns the order of this group, computing it if
        % necessary, as a "vpi" integer
            f = self.orderFun;
            o = f();
        end
        
        function e = elements(self)
        % Returns an Enumerator that enumerates all the elements
        % of this group, using a representation efficient in space
        % if available
            f = self.elementsFun;
            e = f();
        end
        
        function d = decomposition(self)
        % Returns a decomposition of this group, such that every group element
        % can be written as a product of elements from sets
        % (see replab.FiniteGroupDecomposition)
            f = self.decompositionFun;
            d = f();
        end
                
        function g = sample(self)
            g = self.randomBag.sample;
        end
        
        function g = sampleUniformly(self)
        % Returns an element sampled uniformly from this group
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

    properties (Access = protected)
        randomBag_ = []; % Generator for random elements
    end
    
    methods
        
        function R = randomBag(self)
            if isequal(self.randomBag_, [])
                self.randomBag_ = replab.RandomBag(self, self.generators);
            end
            R = self.randomBag_;
        end
        
    end

end
