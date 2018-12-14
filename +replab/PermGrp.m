classdef PermGrp < handle
   
    properties (SetAccess = private)
        n; % domain size
        generators; % m group generators in a m x n matrix, one generator per row
    end
    
    properties (Access = private)
        randomBag_ = [];
        hashedSortedElements_ = [];
        words_ = {};
    end

    methods
        
        function p = randomElement(self)
            if isequal(self.randomBag_, [])
                gens = cell(1, self.nGenerators);
                for i = 1:self.nGenerators
                    gens{i} = self.generator(i);
                end
                self.randomBag_ = replab.prv.RandomBag(gens, 1:self.n, @(x,y) replab.Perm.compose(x, y), @(x) replab.Perm.inverse(x));
            end
            p = self.randomBag_.sample;
        end
        
        function E = allElements(self)
            E = double(self.hashedSortedElements.M');
        end
        
        function o = order(self)
        % Returns a variable precision integer representing the order
        % of this permutation group
            o = vpi(self.hashedSortedElements.nElements);
        end
        
        function p = evaluateWord(self, w)
            p = 1:self.n;
            for i = 1:length(w.indices)
                g = self.generator(w.indices(i));
                e = w.exponents(i);
                we = replab.Perm.pow(g, e);
                p = replab.Perm.compose(p, we);
            end
        end
        
        function N = nGenerators(self)
            N = size(self.generators, 1);
        end
        
        function g = generator(self, i)
            g = self.generators(i, :);
        end
        
        function b = contains(self, permutation)
            b = self.hashedSortedElements.find(permutation) > 0;
        end
        
        function w = factorization(self, permutation)
            ind = self.hashedSortedElements.find(permutation);
            if ind == 0
                error(['Cannot find permutation' num2str(permutation)]);
            end
            w = self.words{ind};
        end
        
    end

    methods (Access = private)
        
        function self = PermGrp(generators)
            self.n = size(generators, 2);
            self.generators = generators;
        end
        
        function b = elementNormalizes(self, el)
            elInv = replab.Perm.inverse(el);
            b = true;
            for i = 1:self.nGenerators
                g = self.generator(i);
                conj = replab.Perm.compose(elInv, replab.Perm.compose(g, el));
                if ~self.contains(conj)
                    b = false;
                    return
                end
            end
        end
        
        function setAllElements(self, hashedSortedElements, words)
            self.hashedSortedElements_ = hashedSortedElements;
            self.words_ = words;
        end
        
        function H = hashedSortedElements(self)
            if isequal(self.hashedSortedElements_, [])
                self.computeAllElements;
            end
            H = self.hashedSortedElements_;
        end
        
        function W = words(self)
            if isequal(self.words_, [])
                self.computeAllElements;
            end
            W = self.words_;
        end
        
        function computeAllElements(self)
            if self.nGenerators == 0
                id = (1:self.n)';
                self.hashedSortedElements_ = replab.prv.HashIntMatrix(id);
                self.words_ = {replab.Word.identity};
            else
                prev = replab.PermGrp(self.generators(1:self.nGenerators-1, :));
                current = prev.closure(self.generators(self.nGenerators, :));
                elements = current.hashedSortedElements;
                [sortedElements indices] = elements.sorted;
                self.hashedSortedElements_ = sortedElements;
                self.words_ = current.words_(indices);
            end
        end
                        
    end
        
    methods (Static)
        
        function G = trivial(n)
            G = replab.PermGrp(zeros(0, n));
        end
        
        function G = symmetric(n)
            if n == 1
                G = replab.PermGrp.trivial(n);
            elseif n == 2
                G = replab.PermGrp.fromGenerators([2 1]);
            else
                G = replab.PermGrp.fromGenerators([2:n 1; 2 1 3:n]);
            end
        end
        
        function G = fromGenerators(permutations)
            n = size(permutations, 2);
            assert(n > 0, 'If creating the trivial group, pass a zeros(0, n) array so that the domain size n is known.');
            m = size(permutations, 1);
            for i = 1:m
                p = permutations(m, :);
                assert(~replab.Perm.isIdentity(p), 'The identity cannot be part of the generators');
            end
            G = replab.PermGrp(permutations);
        end
        
    end
    
    methods % TECHNICAL METHODS
        
        function C = closure(self, el)
            if self.contains(el)
                % trivial closure
                C = self;
                return
            end
            C = replab.PermGrp([self.generators; el]);
            % G is the current group
            % C is the new group
            Gelements = self.hashedSortedElements;
            Gorder = Gelements.nElements;
            Gwords = self.words;
            elWord = replab.Word.generator(C.nGenerators);
            % if <G>^<el> = <G> then <C> = <G> * <el>
            % this follows closely grp.gi in Gap System
            if self.elementNormalizes(el)
                rep = el;
                repWord = elWord;
                shift = 0;
                newElements = zeros(self.n, 0, 'int32');
                newWords = {};
                while ~self.contains(rep)
                    newElements = [newElements zeros(self.n, Gorder)];
                    newWords = horzcat(newWords, cell(1, Gorder));
                    for i = 1:Gorder
                        % we cannot have duplicates here
                        g = Gelements.M(:, i);
                        newElements(:, shift+i) = replab.Perm.compose(g, rep);
                        newWords{shift+i} = Gwords{i} * repWord;
                    end
                    shift = shift + Gorder;
                    rep = replab.Perm.compose(rep, el);
                    repWord = repWord * elWord;
                end
                Celements = Gelements.append(newElements);
                Cwords = horzcat(Gwords, newWords);
                C.setAllElements(Celements, Cwords);
            else
                % otherwise use a Dimino step
                reps = (1:self.n)';
                repWords = {replab.Word.identity};
                Celements = Gelements;
                Cwords = Gwords;
                while size(reps, 2) > 0
                    rep = reps(:, 1);
                    repWord = repWords{1};
                    reps = reps(:, 2:end);
                    repWords = repWords(2:end);
                    for i = 1:C.nGenerators
                        gen = C.generator(i);
                        genWord = replab.Word.generator(i);
                        rg = replab.Perm.compose(rep, gen);
                        rgWord = repWord * genWord;
                        if Celements.find(rg) == 0
                            mem = ismember(rg(:)', Celements.M', 'rows');
                            assert(sum(mem) == 0);
                            newElements = zeros(self.n, Gorder, 'int32');
                            newWords = cell(1, Gorder);
                            for i = 1:Gorder
                                e = Gelements.M(:,i);
                                eWord = Gwords{i};
                                newElements(:,i) = replab.Perm.compose(e, rg);
                                newWords{i} = eWord * rgWord;
                            end
                            reps = [reps rg];
                            repWords = horzcat(repWords, {rgWord});
                            Celements = Celements.append(newElements);
                            Cwords = horzcat(Cwords, newWords);
                        end
                    end
                end
                C.setAllElements(Celements, Cwords);
            end
        end

    end
    
end
