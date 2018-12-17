classdef Chain < handle
% Describes a (generalized) permutation group using a BSGS chain
% constructed using randomized algorithms
    properties
        A; % BSGS action
    end
    
    methods
        
        function cat = G(self)
            cat = self.A.G;
        end
        
        function cat = P(self)
            cat = self.A.P;
        end
        
        % Note: most of those abstract methods are tail-recursive
        % and could be rewritten as loops
        
        function b = isTerm(self)
            b = isa(self, 'replab.bsgs.Term');
        end
        
        function l = length(self)
            l = 0;
            it = self;
            while ~it.isTerm
                l = l + 1;
                it = it.next;
            end
        end
        
        function r = random(self)
            it = self;
            r = 1:self.n;
            while ~it.isTerm
                r = self.A.G.compose(r, it.randomU);
                it = it.next;
            end
        end

        function o = order(self)
            it = self;
            o = vpi(1);
            while ~it.isTerm
                o = o * vpi(it.orbitSize);
                it = it.next;
            end
        end
        
        function remaining = sift(self, el)
            it = self;
            remaining = el;
            while ~it.isTerm
                b = self.A.leftAction(remaining, it.beta);
                i = it.orbitIndex(b);
                if i == 0
                    s = [];
                    return
                else
                    remaining = self.G.compose(it.uInv{i}, remaining);
                end
                it = it.next;
            end
        end
        
        function s = strongGeneratingSet(self)
            s = {};
            it = self;
            while ~it.isTerm
                s = horzcat(s, it.ownSG);
                it = it.next;
            end
        end
        
        function s = strongGeneratingSetWords(self)
            s = {};
            it = self;
            while ~it.isTerm
                s = horzcat(s, it.ownSGWords);
                it = it.next;
            end
        end

        function gd = groupDecomposition(self)
            gd = {};
            it = self;
            while ~it.isTerm
                gd = horzcat(gd, {it.u});
                it = it.next;
            end
        end
        
    end
    
    methods (Static)
      
        function res = siftAndUpdateBaseFrom(A, at, g)
        % Returns an element obtained by sifting g through the chain starting at "at"
        % inserting new base points as required, returns either [],
        % when g can be completely sifted, or
        %
        % {remaining node} where "remaining ~= identity" is the element
        % resulting from an incomplete sift, and node is where
        % this element should be inserted as a strong generator
        %
        % Based on Holt (2005) RANDOMSCHREIER procedure, page 98.
            node = at.next; % we have the chain at -> node
            if node.isTerm
                beta = A.findMovedElement(g);
                if ~isequal(beta, [])
                    newNode = replab.bsgs.Node(A, beta);
                    % insert node so that: at -> newNode -> node
                    at.next = newNode;
                    newNode.next = node;
                    res = {g newNode};
                else
                    res = [];
                end
            else
                b = A.leftAction(g, node.beta);
                i = node.orbitIndex(b);
                if i == 0
                    res = {g node};
                else
                    h = A.G.compose(node.uInv{i}, g);
                    % recursive call
                    res = replab.bsgs.Chain.siftAndUpdateBaseFrom(A, at.next, h);
                end
            end
        end
        
        function addStrongGenerator(start, node, g)
        % Adds a strong generator to a node, and updates the orbits
        % of all its parents
            node.addOwnStrongGenerator(g);
            c = start.next;
            while ~isequal(c, node.next)
                c.update;
                c = c.next;
            end
        end
        
        function b = siftAndAddStrongGenerator(A, start, g)
        % Sifts the given element through the chain and returns whether
        % a new strong generator has been discovered
            res = replab.bsgs.Chain.siftAndUpdateBaseFrom(A, start, g);
            if isequal(res, [])
                b = false;
            else
                g = res{1};
                node = res{2};
                replab.bsgs.Chain.addStrongGenerator(start, node, g);
                b = true;
            end
        end
        
        function chain = randomConstruction(A, randomElement, order, numTests)
        % Constructs a BSGS chain from an oracle that returns random group elements.
        % 
        % If the order of the group is provided, the algorithm is randomized but
        % cannot fail. If the order of the group is not provided (or = -1), then
        % a number of tests numTests is performed so that the probability of failure
        % is 2^-numTests (provided that randomElement returns group elements sampled
        % uniformly at random).
            if nargin < 4
                numTests = 128;
            end
            if nargin < 3 || isequal(order, [])
                hasOrder = false;
            else
                hasOrder = true;
            end
            start = replab.bsgs.Start.emptyChain(A);
            numSifted = 0;
            while (~hasOrder && numSifted <= numTests) || (hasOrder && start.next.order < order)
                s = randomElement();
                b = replab.bsgs.Chain.siftAndAddStrongGenerator(A, start, s);
                if b
                    numSifted = 0;
                else
                    numSifted = numSifted + 1;
                end
            end
            chain = start.next;
        end
        
        function chain = fromGenerators(A, generators, order, numTests)
        % Constructs a BSGS chain from a list of generators, given a row vectors
        % in a matrix.
        %
        % If the order of the group is provided, the algorithm is randomized but
        % cannot fail.
        %
        % Random elements are produced using the product replacement algorithm;
        % that algorithm can be unsatisfactory when the group is a direct product
        % of a large number of copies of the same finite simple group 
        % (see Holt et al. Handbook of Computational Group Theory (2005), p. 69).
            if nargin < 4
                numTests = 128;
            end
            if nargin < 3 || isequal(order, [])
                order = [];
            end
            bag = replab.prv.RandomBag(A.G, generators);
            chain = replab.bsgs.Chain.randomConstruction(A, @() bag.sample, order, numTests);
        end
        
    end
    
end
