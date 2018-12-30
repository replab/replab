classdef Chain < handle
% Describes a (generalized) permutation group using a BSGS chain
% constructed using randomized algorithms
    properties (Access = protected)
        A; % BSGS action
    end
    
    methods
        
        function check(self, group)
            it = self;
            A = self.A;
            P = A.P;
            G = A.G;
            while ~it.isTerm
                beta = it.beta;
                it1 = it.next;
                while ~it1.isTerm
                    for i = 1:length(it1.orbit)
                        P.assertEqv(A.leftAction(it1.u{i}, beta), beta);
                    end
                    it1 = it1.next;
                end
                for i = 1:length(it.orbit)
                    u = it.u{i};
                    uInv = it.uInv{i};
                    uInvW = it.uInvW{i};
                    G.assertEqv(uInv, group.evaluateWord(uInvW));
                    G.assertEqv(G.compose(u, uInv), G.identity);
                    b = it.orbit{i};
                    P.assertEqv(A.leftAction(u, beta), b);
                    P.assertEqv(A.leftAction(uInv, b), beta);
                end
                it = it.next;
            end
        end
        
        function b = areWordsCompleted(self)
            it = self;
            b = true;
            while ~it.isTerm
                if it.emptyWords > 0
                    b = false;
                    return
                end
                it = it.next;
            end
        end
        
        function [r rw] = wordsStep(self, t, tw)
            assert(~self.isTerm);
            G = self.G;
            A = self.A;
            b = A.leftAction(t, self.beta);
            ind = self.orbitIndex(b);
            assert(ind ~= 0);
            uInvW = self.uInvW{ind};
            if ~isequal(uInvW, [])
                r = G.compose(self.uInv{ind}, t);
                rw = uInvW * tw;
                if tw.length < uInvW.length
                    twInv = inv(tw);
                    tInv = G.inverse(t);
                    self.uInvW{ind} = twInv;
                    self.uInv{ind} = tInv;
                    self.u{ind} = G.inverse(tInv);
                    self.newWord(ind) = true;
                    self.wordsStep(tInv, twInv);
                end
            else
                twInv = inv(tw);
                tInv = G.inverse(t);
                r = G.identity;
                rw = replab.Word.identity;
                self.uInvW{ind} = twInv;
                self.uInv{ind} = tInv;
                self.u{ind} = G.inverse(tInv);
                self.newWord(ind) = true;
                self.emptyWords = self.emptyWords - 1;
                self.wordsStep(tInv, twInv);
            end
        end
        
        function [t tw] = wordsRound(self, l, t, tw)
            G = self.G;
            it = self;
            while ~it.isTerm && ~G.isIdentity(t) && tw.length < l
                [t tw] = it.wordsStep(t, tw);
                it = it.next;
            end
        end
        
        function wordsImprove(self, l)
            G = self.G;
            it = self;
            while ~it.isTerm
                nonEmpty = [];
                for i = 1:it.orbitSize
                    if ~isequal(it.uInvW{i}, [])
                        nonEmpty = [nonEmpty i];
                    end
                end
                for i = nonEmpty
                    for j = nonEmpty
                        if it.newWord(i) || it.newWord(j)
                            t = G.compose(it.uInv{i}, it.uInv{j});
                            tw = it.uInvW{i} * it.uInvW{j};
                            it.wordsRound(l, t, tw);
                        end
                    end
                end
                it.newWord(:) = false;
                it = it.next;
            end
        end
        
        function wordsQuick(self, group, maxLength, s, l)
            count = 1;
            nGens = length(group.generators);
            G = self.G;
            for len = 1:maxLength
                numWords = nGens^len;
                indices = zeros(1, len);
                for count = 1:numWords
                    t = G.identity;
                    rem = count - 1;
                    for i = 1:len
                        g = mod(rem, nGens);
                        rem = (rem - g)/nGens;
                        indices(i) = g + 1;
                        t = G.compose(t, group.generators{g+1});
                    end
                    tw = replab.Word.fromIndicesAndExponents(indices, ones(1, len));
                    G.assertEqv(t, group.evaluateWord(tw));
                    self.wordsRound(l, t, tw);
                    if mod(count, s) == 0
                        self.wordsImprove(l);
                        l = l * 5/4;
                    end
                    if self.areWordsCompleted
                        return
                    end
                end
            end
        end
        
        function wordsComplete(self)
        % Completes the word decomposition in the BSGS chain
            if ~self.isTerm
                self.next.wordsComplete;
                A = self.A;
                G = self.G;
                while self.emptyWords > 0
                    self
                    for i = 1:self.orbitSize
                        if isempty(self.uInvW{i})
                            b = self.orbit{i};
                            it = self;
                            solved = false;
                            while ~solved && ~it.isTerm
                                for j = 1:it.orbitSize
                                    tw = it.uInvW{j};
                                    if ~isempty(tw)
                                        t = it.uInv{j};
                                        newB = A.leftAction(t, b);
                                        ind = self.orbitIndex(newB);
                                        assert(ind > 0);
                                        if ~isempty(self.uInvW{ind})
                                            uInvW = self.uInvW{ind} * tw;
                                            uInv = G.compose(self.uInv{ind}, t);
                                            self.uInvW{i} = uInvW;
                                            self.uInv{i} = uInv;
                                            self.u{i} = G.inverse(uInv);
                                            self.emptyWords = self.emptyWords - 1;
                                            solved = true;
                                            break
                                        end
                                    end
                                end
                                it = it.next;
                            end
                        end
                    end
                end
            end
        end
        
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
        
        function s = orbitSizes(self)
            s = zeros(1, self.baseLength);
            it = self;
            i = 1;
            while ~it.isTerm
                s(i) = it.orbitSize;
                i = i + 1;
                it = it.next;
            end
        end
        
        function l = baseLength(self)
            l = 0;
            it = self;
            while ~it.isTerm
                l = l + 1;
                it = it.next;
            end
        end
        
        function b = base(self)
            b = [];
            it = self;
            while ~it.isTerm
                b(end+1) = it.beta;
                it = it.next;
            end
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
            r = self.G.identity;
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
        
        function [word remaining] = factor(self, el)
            word = [];
            [remaining indices] = self.sift(el);
            if ~isempty(indices)
                word = self.wordFromIndices(indices);
            end
        end
        
        function el = elementFromIndices(self, indices)
            G = self.G;
            i = 1;
            it = self;
            el = G.identity;
            while ~it.isTerm
                el = G.compose(el, it.u{indices(i)});
                it = it.next;
                i = i + 1;
            end
        end
        
        function word = wordFromIndices(self, indices)
            i = 1;
            it = self;
            word = replab.Word.identity;
            while ~it.isTerm
                word = word / it.uInvW{indices(i)};
                it = it.next;
                i = i + 1;
            end
        end
        
        function factors = indexFactors(self)
            s = self.orbitSizes;
            factors = fliplr(cumprod(vpi(fliplr(s))));
            factors = [factors(2:end) vpi(1)];
        end
        
        function indices = indicesFromIndex(self, index)
            L = self.baseLength;
            indices = zeros(1, L);
            f = self.orbitSizes;
            ind = index - 1;
            for i = L:-1:1
                r = mod(ind, f(i));
                ind = (ind - r)/f(i);
                indices(i) = double(r) + 1;
            end
        end
        
        function index = indexFromIndices(self, indices)
            index = vpi(0);
            L = self.baseLength;
            f = self.orbitSizes;
            for i = 1:L
                index = index * f(i);
                index = index + vpi(indices(i) - 1);
            end
            index = index + vpi(1);
        end

        function [remaining indices] = sift(self, el)
            it = self;
            remaining = el;
            indices = zeros(1, self.baseLength);
            i = 1;
            while ~it.isTerm
                b = self.A.leftAction(remaining, it.beta);
                ind = it.orbitIndex(b);
                if ind == 0
                    remaining = [];
                    return
                else
                    indices(i) = ind;
                    remaining = self.G.compose(it.uInv{ind}, remaining);
                end
                i = i + 1;
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
        
        function chain = forBSGSGroup(group)
        % Constructs a BSGS chain for a BSGS group.
        %
        % If the order of the group is provided, the algorithm is randomized but
        % cannot fail.
        %
        % Random elements are produced using the product replacement algorithm;
        % that algorithm can be unsatisfactory when the group is a direct product
        % of a large number of copies of the same finite simple group 
        % (see Holt et al. Handbook of Computational Group Theory (2005), p. 69).
            numTests = ceil(-log2(replab.Settings.bsgsFailureProbability));
            if group.knownOrder
                order = group.order;
            else
                order = [];
            end
            bag = replab.RandomBag(group, group.generators);
            chain = replab.bsgs.Chain.randomConstruction(group.action, @() bag.sample, order, numTests);
        end
        
    end
    
end
