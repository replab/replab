% Generator of random elements from generators of a group.
% 
% A random bag is a set of random group elements that always generates
% the group; random elements are provided by multiplying elements of the
% bag and returning one element of the product which is removed from the bag.
% 
% Straight-forward implementation of PRINITIALIZE and PRRANDOM of 
% section 3.2.2, pp. 70-71 of Holt 2005 (Handbook of Computational Group Theory)
%
% Is generic in the group element type, using the replab.cat framework.
classdef RandomBag < handle
    
    properties
        G % Group definition
        x0 % Last generated sample
        x % 1 x r cell array representing the contents of the bag
    end

    methods
        
        function res = sample(self)
            r = length(self.x);
            s = randi(r);
            t = randi(r);
            while t == s
                t = randi(r);
            end
            if randi(2) == 2
                if randi(2) == 2 % e = 1
                    self.x{s} = self.G.compose(self.x{s}, self.x{t});
                else
                    self.x{s} = self.G.composeWithInverse(self.x{s}, self.x{t});
                end
                self.x0 = self.G.compose(self.x0, self.x{s});
            else
                if randi(2) == 2 % e = 1
                    self.x{s} = self.G.compose(self.x{t}, self.x{s});
                else
                    tinv = self.G.inverse(self.x{t});
                    self.x{s} = self.G.compose(tinv, self.x{s});
                end
                self.x0 = self.G.compose(self.x{s}, self.x0);
            end
            res = self.x0;
        end
        
        function self = RandomBag(G, generators, r, n)
        % Constructs a random bag from the given generators, given
        % as a 1 x k cell array of group elements, where k >= 0.
        %
        % self.cat must be an instance of replab.prv.Group
        %
        % r is the number of elements in the bag (optional, self.catault: max(k, 10))
        % n is the number of shuffles done during initialization (optional, self.catault: 50)
            if nargin < 4
                r = -1;
            end
            if nargin < 3
                n = 50;
            end
            k = length(generators); % number of generators
            if r < k || r < 10
                r = max(k, 10);
            end
            x = cell(1, r);
            if k == 0
                % cater for the special case when generators are empty
                for i = 1:r
                    x{i} = G.identity;
                end
            else
                g = 1;
                for i = 1:r
                    x{i} = generators{g};
                    g = g + 1;
                    if g > k
                        g = 1;
                    end
                end
            end
            self.G = G;
            self.x0 = G.identity; % initially, the identity element
            self.x = x;
            for i = 1:n
                self.sample; % perform initial shuffles
            end
        end
   
    end
    
end
