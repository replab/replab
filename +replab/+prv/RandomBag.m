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
        withWords
        word0
        words
    end

    methods
        
        function [res resWord] = sample(self)
            r = length(self.x);
            s = randi(r);
            t = randi(r);
            while t == s
                t = randi(r);
            end
            if randi(2) == 2
                if randi(2) == 2 % e = 1
                    self.x{s} = self.G.compose(self.x{s}, self.x{t});
                    if self.withWords
                        self.words{s} = self.words{s} * self.words{t};
                    end
                else
                    self.x{s} = self.G.composeWithInverse(self.x{s}, self.x{t});
                    if self.withWords
                        self.words{s} = self.words{s} * inv(self.words{t});
                    end
                end
                self.x0 = self.G.compose(self.x0, self.x{s});
                if self.withWords
                    self.word0 = self.word0 * self.words{s};
                end
            else
                if randi(2) == 2 % e = 1
                    self.x{s} = self.G.compose(self.x{t}, self.x{s});
                    if self.withWords
                        self.words{s} = self.words{t} * self.words{s};
                    end
                else
                    tinv = self.G.inverse(self.x{t});
                    self.x{s} = self.G.compose(tinv, self.x{s});
                    if self.withWords
                        self.words{s} = inv(self.words{t}) * self.words{s};
                    end
                end
                self.x0 = self.G.compose(self.x{s}, self.x0);
                if self.withWords
                    self.word0 = self.words{s} * self.word0;
                end                
            end
            res = self.x0;
            if self.withWords
                resWord = self.word0;
            end
        end
        
        function self = RandomBag(G, generators, withWords, r, n)
        % Constructs a random bag from the given generators, given
        % as a 1 x k cell array of group elements, where k >= 0.
        %
        % self.cat must be an instance of replab.prv.Group
        %
        % r is the number of elements in the bag (optional, self.catault: max(k, 10))
        % n is the number of shuffles done during initialization (optional, self.catault: 50)
            if nargin < 5
                r = -1;
            end
            if nargin < 4
                n = 50;
            end
            if nargin < 3
                withWords = false;
            end
            k = length(generators); % number of generators
            if r < k || r < 10
                r = max(k, 10);
            end
            x = cell(1, r);
            if withWords
                words = cell(1, r);
            else
                words = [];
            end
            if k == 0
                % cater for the special case when generators are empty
                for i = 1:r
                    x{i} = G.identity;
                end
                if withWords
                    words{i} = replab.Word.identity;
                end
            else
                g = 1;
                for i = 1:r
                    x{i} = generators{g};
                    if withWords
                        words{i} = replab.Word.generator(g);
                    end
                    g = g + 1;
                    if g > k
                        g = 1;
                    end
                end
            end
            self.G = G;
            self.x0 = G.identity; % initially, the identity element
            if withWords
                self.word0 = replab.Word.identity;
            else
                self.word0 = [];
            end
            self.x = x;
            self.withWords = withWords;
            self.words = words;
            for i = 1:n
                self.sample; % perform initial shuffles
            end
        end
   
    end
    
end
