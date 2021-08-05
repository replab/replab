classdef PermutationGroup < replab.FiniteGroup
% A base class for all permutation groups

    properties (SetAccess = protected)
        domainSize % integer: The integer ``d``, as this group acts on ``{1, ..., d}``
    end

    methods % Constructor

        function self = PermutationGroup(domainSize, generators, varargin)
        % Constructs a permutation group
        %
        % Args:
        %   domainSize (integer): Size of the domain
        %   generators (cell(1,\*) of permutation): Group generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi, optional): Order of the group
        %   type (`+replab.PermutationGroup`, optional): Type of this group if known,
        %                                                or ``'self'`` if this group is its own type
        %   relators (cell(1,\*) of charstring): Relators
        %   chain (`+replab.+bsgs.Chain`): BSGS chain describing the group
            identity = 1:domainSize;
            for i = 1:length(generators)
                assert(~all(generators{i} == identity), 'Generators cannot contain the identity');
            end
            args = struct('type', [], 'chain', [], 'relators', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            if ~isempty(args.type)
                type = args.type;
            else
                type = replab.S(domainSize);
            end
            self@replab.FiniteGroup(identity, generators, type, restArgs{:});
            self.domainSize = domainSize;
            if ~isempty(args.chain)
                assert(~args.chain.isMutable);
                self.cache('chain', args.chain, 'ignore');
            end
        end

    end

    methods (Static)

        function pg = fromChain(chain, type)
            if nargin < 2
                type = [];
            end
            pg = replab.PermutationGroup(chain.n, chain.strongGenerators, 'order', chain.order, 'type', type, 'chain', chain);
        end

    end

    methods (Access = protected)

        function f = computeFactorization(self)
            f = replab.mrp.Factorization.make(self);
        end

        function c = computeLexChain(self)
            c = self.chain;
            if ~c.isLex
                c = c.mutableCopy;
                c.baseChange(1:self.domainSize, true);
                c.makeImmutable;
            end
        end

        function chain = computeChain(self)
            chain = self.cachedOrEmpty('partialChain');
            if isempty(chain)
                chain = replab.bsgs.Chain.make(self.domainSize, self.generators, [], self.cachedOrEmpty('order'));
            else
                if chain.isMutable
                    if self.inCache('order')
                        chain.randomizedSchreierSims(self.cachedOrEmpty('order'));
                    else
                        chain.deterministicSchreierSims;
                    end
                    chain.makeImmutable;
                    self.cache('partialChain', chain, 'overwrite');
                end
            end
        end

        function c = computePartialChain(self)
            if self.inCache('chain')
                c = self.chain;
            else
                c = replab.bsgs.Chain.makeBoundedOrder(self.domainSize, self.generators, replab.globals.fastChainOrder);
            end
        end

    end

    methods (Access = protected) % Implementations

        function o = computeOrder(self)
            o = self.chain.order;
        end

        function E = computeElements(self)
            basis = replab.util.MixedRadix(self.lexChain.orbitSizes, true, false);
            atFun = @(ind) self.lexChain.elementFromIndices(basis.ind2sub(ind));
            findFun = @(el) basis.sub2ind(self.lexChain.indicesFromElement(el));
            E = replab.IndexedFamily.lambda(self.order, atFun, findFun);
        end

        function classes = computeConjugacyClasses(self)
            if self.order < 100000
                C = replab.perm.conjugacyClassesByOrbits(self);
                n = length(C);
                classes = cell(1, n);
                for i = 1:n
                    cl = sortrows(C{i}');
                    classes{i} = replab.ConjugacyClass(self, cl(1,:));
                end
            else
                classes = replab.perm.conjugacyClassesByRandomSearch(self);
            end
            reps = zeros(length(classes), self.domainSize);
            for i = 1:length(classes)
                reps(i,:) = classes{i}.representative;
            end
            [~, I] = sortrows(reps);
            classes = replab.ConjugacyClasses(self, classes(I)); % sort by minimal representative
        end

        function res = computeIsCyclic(self)
            if self.nGenerators <= 1
                res = true;
            elseif ~self.isCommutative
                res = false;
            else
                pds = factor(self.order);
                if length(pds) == 1 % group of prime order is cyclic
                    res = true;
                    return
                end
                assert(all(pds <= 2^53-1)); % to be sure (otherwise can a BSGS have been computed?)
                pds = unique(double(pds));
                pds = double(pds);
                for p = pds
                    newGens = cellfun(@(g) self.composeN(g, p), self.generators, 'uniform', 0);
                    newGens = newGens(~cellfun(@(g) self.isIdentity(g), newGens));
                    if self.subgroup(newGens).order*p ~= self.order
                        res = false;
                        return
                    end
                end
                res = true;
            end
        end

        function res = computeIsSimple(self)
            if self.isTrivial
                res = false;
                return
            end
            C = self.conjugacyClasses;
            for i = 1:C.nClasses
                c = C.classes{i}.representative;
                if ~self.isIdentity(c)
                    if self.normalClosure(self.subgroup({c})) ~= self
                        res = false;
                        return
                    end
                end
            end
            res = true; % all conjugacy classes generate the full group
        end

        function sub = computeDerivedSubgroup(self)
            nG = self.nGenerators;
            n = self.domainSize;
            chain = replab.bsgs.Chain(n);
            generators = {};
            for i = 1:nG
                gi = self.generator(i);
                for j = 1:nG
                    gj = self.generator(j);
                    cm = self.composeWithInverse(self.compose(gi, gj), self.compose(gj, gi));
                    if chain.stripAndAddStrongGenerator(cm)
                        generators{1, end+1} = cm;
                        chain.randomizedSchreierSims([]);
                    end
                end
            end
            % compute the normal closure
            toCheck = generators;
            while ~isempty(toCheck)
                h = toCheck{end};
                toCheck = toCheck(1:end-1);
                for i = 1:nG
                    gi = self.generator(i);
                    cm = self.leftConjugate(gi, h);
                    if chain.stripAndAddStrongGenerator(cm)
                        generators{1, end+1} = cm;
                        toCheck{1, end+1} = cm;
                        chain.randomizedSchreierSims([]);
                    end
                end
            end
            chain.makeImmutable;
            sub = replab.PermutationGroup(self.domainSize, generators, 'order', chain.order, 'type', self.type, 'chain', chain);
        end

        function s = computeSetProduct(self)
            c = self.chain;
            k = c.length;
            T = cell(1, k);
            for i = 1:k
                Ui = c.U{i};
                m = size(Ui, 2);
                Ti = cell(1, m);
                for j = 1:m
                    Ti{j} = Ui(:,j)';
                end
                T{i} = Ti;
            end
            s = replab.SetProduct(self, T, true);
        end

    end

    methods (Access = protected)

        function G = computeNiceGroup(self)
            G = self;
        end

        function m = computeNiceMorphism(self)
            m = replab.FiniteIsomorphism.identity(self);
        end

    end

    methods % Group internal description

        function c = chain(self)
        % Returns the stabilizer chain corresponding to this permutation group.
        %
        % Returns:
        %   `+replab.+bsgs.Chain`: Stabilizer chain
            c = self.cached('chain', @() self.computeChain);
        end

        function f = factorization(self)
        % Returns an object able to compute factorizations in the group generators
        %
        % Returns:
        %   `+replab.+mrp.Factorization`: A factorization instance
            f = self.cached('factorization', @() self.computeFactorization);
        end

        function c = lexChain(self)
        % Returns the reduced stabilizer chain corresponding to this permutation group in lexicographic order
        %
        % No base points are redundant.
        %
        % Returns:
        %   `+replab.+bsgs.Chain`: Stabilizer chain
            c = self.cached('lexChain', @() self.computeLexChain);
        end

        function c = partialChain(self)
        % Returns the stabilizer chain corresponding to this permutation group if it can be computed quickly
            c = self.cached('partialChain', @() self.computePartialChain);
        end

    end

    methods (Access = protected)

        % Morphisms

        function m = morphismByImages_(self, target, preimages, images, imageElementFun)
            assert(length(preimages) == length(images));
            if isa(target, 'replab.PermutationGroup')
                m = replab.mrp.PermToPerm(self, target, preimages, images, imageElementFun);
            elseif isa(target, 'replab.FiniteGroup')
                m = replab.mrp.PermToFiniteGroup(self, target, preimages, images, imageElementFun);
            else
                m = replab.mrp.PermToGroup(self, target, preimages, images, imageElementFun);
            end
        end

    end

    methods % Implementations

        % replab.Str

        function s = headerStr(self)
            if self.inCache('order')
                s = sprintf('Permutation group acting on %d elements of order %s', self.domainSize, strtrim(num2str(self.order)));
            else
                s = sprintf('Permutation group acting on %d elements', self.domainSize);
            end
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.FiniteGroup(self);
            names{1,end+1} = 'domainSize';
        end

        % replab.Obj

        function l = laws(self)
            l = replab.laws.PermutationGroupLaws(self);
        end

        % replab.Domain

        function b = eqv(self, x, y)
            b = isequal(x, y);
        end

        function s = sample(self)
            s = self.chain.sample;
        end

        % replab.Monoid

        function z = compose(self, x, y)
            z = x(y);
        end

        % replab.Group

        function y = inverse(self, x)
            n = length(x);
            y = zeros(1, n);
            y(x) = 1:n;
        end

        function z = composeWithInverse(self, x, y)
            z = zeros(1, length(x));
            z(y) = x;
        end

        function z = leftConjugate(self, x, y)
            z = zeros(1, length(x));
            % x y xInv
            z(x) = x(y);
        end

        % FiniteGroup

        function res = withGeneratorNames(self, newNames)
        % Returns a modified copy of this finite group with renamed generators
        %
        % Args:
        %   newNames (cell(1,\*) of charstring): New generator names
        %
        % Returns:
        %   `.FiniteGroup`: Updated copy
            if isequal(self.generatorNames, newNames)
                res = self;
                return
            end
            args = cell(1, 0);
            if self.inCache('order')
                args = horzcat(args, {'order', self.order});
            end
            if self.inCache('chain')
                args = horzcat(args, {'chain', self.chain});
            end
            res = replab.PermutationGroup(self.domainSize, self.generators, 'type', self.type, 'generatorNames', newNames, args{:});
        end

        function A = abstractGroup(self, generatorNames)
            if nargin < 2 || isempty(generatorNames)
                generatorNames = self.generatorNames;
            end
            A = replab.AbstractGroup(generatorNames, self.relators(generatorNames), 'permutationGenerators', self.generators);
        end

        function m = abstractMorphism(self, generatorNames)
            if nargin < 2 || isempty(generatorNames)
                generatorNames = self.generatorNames;
            end
            m = self.abstractGroup(generatorNames).niceMorphism.inverse;
        end

    end

    methods % Group elements

        function b = contains(self, g)
        % Tests whether this group contains the given element
        %
        % Args:
        %   g (element of `.type`): Element to test membership of
        %
        % Returns:
        %   logical: True if this group contains ``g`` and false otherwise
            b = self.chain.contains(g);
        end

        function o = elementOrder(self, p)
            o = replab.Permutation.order(p);
        end


        function g = imageLetters(self, letters)
            g = self.identity;
            L = length(letters);
            for i = 1:L
                l = letters(i);
                if l > 0
                    g = g(self.generator(l));
                else
                    g(self.generator(-l)) = g;
                end
            end
        end

        function l = factorizeLetters(self, element)
            l = self.factorization.factorize(element);
        end

    end

    methods % Construction of groups

        function res = closure(self, rhs)
            if isa(rhs, 'replab.PermutationGroup')
                % if one group contains the other
                if self.isSubgroupOf(rhs)
                    res = rhs;
                    return
                end
                if rhs.isSubgroupOf(self)
                    res = self;
                    return
                end
                % otherwise do the computation
                c = self.chain.mutableCopy;
                for i = 1:rhs.nGenerators
                    if c.stripAndAddStrongGenerator(rhs.generator(i))
                        c.randomizedSchreierSims([]);
                    end
                end
            else
                % if the group already contains the element
                if self.contains(rhs)
                    res = self;
                    return
                end
                % otherwise do the computation
                c = self.chain.mutableCopy;
                if c.stripAndAddStrongGenerator(rhs)
                    c.randomizedSchreierSims([]);
                end
            end
            c.makeImmutable;
            res = replab.PermutationGroup.fromChain(c, self.type);
        end

        function nc = normalClosure(self, rhs)
            if isa(rhs, 'replab.FiniteGroup')
                toCheck = rhs.generators;
            else
                toCheck = {rhs};
            end
            chain = replab.bsgs.Chain(self.domainSize);
            generators = {};
            while ~isempty(toCheck)
                rhsg = toCheck{end};
                toCheck = toCheck(1:end-1);
                for i = 1:self.nGenerators
                    gi = self.generator(i);
                    cm = self.leftConjugate(gi, rhsg);
                    if chain.stripAndAddStrongGenerator(cm)
                        generators{1, end+1} = cm;
                        toCheck{1, end+1} = cm;
                        chain.randomizedSchreierSims([]);
                    end
                end
            end
            chain.makeImmutable;
            nc = replab.PermutationGroup(self.domainSize, generators, 'order', chain.order, 'type', self.type, 'chain', chain);
        end

        % Subgroups

        function sub = subgroupWithGenerators(self, generators, order)
        % Constructs a permutation subgroup from its generators
        %
        % Args:
        %   generators (cell array): List of generators given as a permutations in a row cell array
        %   order (vpi or ``[]``, optional): Argument specifying the group order, if given speeds up computations
        %
        % Returns:
        %   +replab.PermutationGroup: The constructed permutation subgroup
            if nargin < 3
                order = [];
            end
            if length(generators) == self.nGenerators && all(arrayfun(@(i) self.eqv(self.generator(i), generators{i}), 1:length(generators)))
                sub = self;
            else
                sub = replab.PermutationGroup(self.domainSize, generators, 'order', order, 'type', self.type);
            end
        end

        function c = centralizer(self, other)
            if ~isa(other, 'replab.PermutationGroup')
                other = self.subgroup({other});
            end
            c = replab.bsgs.Centralizer(self, other).subgroup;
        end

        function res = intersection(self, other)
            assert(self.hasSameTypeAs(other));
            if self.order > other.order
                res = replab.bsgs.Intersection(other, self).subgroup;
            else
                res = replab.bsgs.Intersection(self, other).subgroup;
            end
        end

        % Cosets

        function c = findLeftConjugations(self, s, t, sCentralizer, tCentralizer)
            if nargin < 4 || isempty(sCentralizer)
                sCentralizer = self.centralizer(s);
            end
            if nargin < 5 || isempty(tCentralizer)
                tCentralizer = [];
            end
            b = replab.bsgs.LeftConjugation(self, s, t, sCentralizer, tCentralizer).find;
            % note that we have
            % sCentralizer == tCentralizer.leftConjugateGroup(self.inverse(b))
            % tCentralizer == sCentralizer.leftConjugateGroup(b)
            if isempty(b)
                c = [];
            else
                c = sCentralizer.leftCoset(b, self);
            end
        end

        % Relations to other groups

        function res = hasSameTypeAs(self, rhs)
            res = isa(rhs, 'replab.PermutationGroup') && (self.type.domainSize == rhs.type.domainSize);
        end

        % Representations

        function rep = regularRep(self)
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
            rep = self.permutationRep(o, 'preimages', self.generators, 'images', perms);
        end

    end

    methods % Methods specific to permutation groups

        function G = generatorsAsMatrix(self)
        % Returns the generators of this group concatenated in a matrix
        %
        % Returns:
        %   integer(\*,\*): Matrix whose rows are the group generators
            G = zeros(self.nGenerators, self.domainSize);
            for i = 1:self.nGenerators
                G(i,:) = self.generator(i);
            end
        end

        function sub = unorderedPartitionStabilizer(self, partition)
        % Computes the subgroup that leaves the given unordered partition invariant
        %
        % Example:
        %    >>> S4 = replab.S(4);
        %    >>> G = S4.unorderedPartitionStabilizer(replab.Partition.fromVector([1 1 2 2]));
        %    >>> G == S4.subgroup({[2 1 3 4] [3 4 1 2]})
        %        1
        %
        % Args:
        %   partition (`.Partition`): Unordered partition to leave invariant
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup that stabilizes the unordered partition
            sub = replab.bsgs.UnorderedPartitionStabilizer(self, partition).subgroup;
        end

        function sub = orderedPartitionStabilizer(self, partition)
        % Computes the subgroup that leaves the given ordered partition invariant
        %
        % The subgroup maps every block of the partition to itself under permutation of the
        % domain elements.
        %
        %
        % Example:
        %    >>> S4 = replab.S(4);
        %    >>> G = S4.orderedPartitionStabilizer(replab.Partition.fromVector([1 1 2 2]));
        %    >>> G == S4.subgroup({[2 1 3 4] [1 2 4 3]})
        %        1
        %
        % Args:
        %   partition (`.Partition`): Ordered partition to leave invariant
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup that stabilizes the ordered partition
            sub = replab.bsgs.OrderedPartitionStabilizer(self, partition).subgroup;
        end

        function [sMinLex, P] = vectorFindLexMinimal(self, s, sStabilizer)
        % Finds the lexicographic minimal vector under permutation of its coefficients by this group
        %
        % Example:
        %   >>> n = 10;
        %   >>> s = randi([-3 3], 1, n);
        %   >>> G = replab.S(n);
        %   >>> [sml, P] = G.vectorFindLexMinimal(s);
        %   >>> all(sml == sort(s)) % lexmin using the symmetric group is simply a sort
        %       1
        %   >>> all(sml(P.representative) == s)
        %       1
        %
        % Args:
        %   s (double(1,domainSize)): Vector to permute
        %   sStabilizer (`.PermutationGroup` or ``[]``, optional): Stabilzier of ``s``
        %
        % Returns
        % -------
        %   sMinLex: double(1,domainSize):
        %     Minimal lexicographic representative of ``s`` under permutation by this group
        %   P: `.LeftCoset`
        %     Set of permutations ``p`` such that ``sMinLex == s(inverse(p))`` or ``s == sMinLex(p)``
            if nargin < 3 || isempty(sStabilizer)
                sStabilizer = self.vectorStabilizer(s);
            end
            lm = replab.bsgs.LexMin(self, s, sStabilizer);
            [sMinLex, p] = lm.search;
            P = sStabilizer.leftCoset(p, self);
        end

        function P = vectorFindPermutationsTo(self, s, t, sStabilizer, tStabilizer)
        % Finds the permutations that send a vector to another vector
        %
        % We return the set of ``p`` such that ``t == s(inverse(p))`` or ``s == t(p)``.
        %
        % We use this notation as the left action of ``p`` on a list ``s`` is given by ``s(inverse(p))``.
        %
        % Args:
        %   s (double(1,domainSize)): Source vector
        %   t (double(1,domainSize)): Target vector
        %   sStabilizer (`.PermutationGroup` or ``[]``, optional): Stabilizer of ``s``
        %   tStabilizer (`.PermutationGroup` or ``[]``, optional): Stabilizer of ``t``
        %
        % Returns:
        %   `+replab.LeftCoset`: The set of permutations ``{p}`` such that ``t == s(inverse(p))``; or ``[]`` if no element found
            if nargin < 4 || isequal(sStabilizer, [])
                sStabilizer = self.vectorStabilizer(s);
            end
            if nargin < 5 || isequal(tStabilizer, [])
                tStabilizer = self.vectorStabilizer(t);
            end
            p = replab.bsgs.VectorPermutationTo(self, s, t, sStabilizer, tStabilizer).find;
            if ~isempty(p)
                P = sStabilizer.leftCoset(p, self);
            else
                P = [];
            end
        end

        function P = matrixFindPermutationsTo(self, S, T, SStabilizer, TStabilizer)
        % Finds the permutations that send a matrix to another matrix
        %
        % We return the set of ``p`` such that ``T == S(inverse(p),inverse(p))`` or ``S == T(p,p)``.
        %
        % We use this notation as the left action of ``p`` on a matrix ``S`` is given by ``S(inverse(p),inverse(p))``.
        %
        % Args:
        %   S (double(1,domainSize)): Source matrix
        %   T (double(1,domainSize)): Target matrix
        %   SStabilizer (`.PermutationGroup` or ``[]``, optional): Stabilizer of ``S``
        %   TStabilizer (`.PermutationGroup` or ``[]``, optional): Stabilizer of ``T``
        %
        % Returns:
        %   `+replab.LeftCoset`: The set of permutations ``{p}`` such that ``T == S(inverse(p),inverse(p))``; or ``[]`` if no element found
            if nargin < 4 || isequal(SStabilizer, [])
                SStabilizer = self.matrixStabilizer(S);
            end
            if nargin < 5 || isequal(TStabilizer, [])
                TStabilizer = self.matrixStabilizer(T);
            end
            p = replab.bsgs.MatrixPermutationTo(self, S, T, SStabilizer, TStabilizer).find;
            if ~isempty(p)
                P = SStabilizer.leftCoset(p, self);
            else
                P = [];
            end
        end

        function sub = vectorStabilizer(self, vector)
        % Returns the permutation subgroup that leaves a given vector invariant
        %
        % Args:
        %   vector (double(1,domainSize)): Vector to stabilize under permutation
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup of this group leaving ``vector`` invariant
            partition = replab.Partition.fromVector(vector);
            sub = self.orderedPartitionStabilizer(partition);
        end

        function sub = matrixStabilizer(self, matrix)
        % Returns the permutation subgroup that leaves a given matrix invariant by joint permutation of its rows and columns
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> A = [1 1 0 1; 1 1 1 0; 0 1 1 1; 1 0 1 1]; % adjacency matrix of the square
        %   >>> G = S4.matrixStabilizer(A);
        %   >>> G.order == 8 % is the dihedral group of order 8
        %       1
        %
        % Args:
        %   matrix (double(domainSize,domainSize)): Matrix to stabilize under permutation
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup of this group such that every element ``g`` satisfies ``matrix(g,g) == matrix``
            sub = replab.bsgs.MatrixAutomorphism(self, matrix).subgroup;
        end

        function s = setwiseStabilizer(self, set)
        % Returns the subgroup that stabilizes the given set as a set
        %
        % i.e. for this group ``G``, it returns ``H = {g \in G : g(set) = set}``
        %
        % Example:
        %   >>> G = replab.S(4).subgroup({[3 1 2 4] [1 4 2 3]});
        %   >>> H = G.setwiseStabilizer([1 2]);
        %   >>> H == replab.S(4).subgroup({[2 1 4 3]})
        %       1
        %
        % Example:
        %   >>> G = replab.S(10).subgroup({[1,3,2,10,9,8,6,5,7,4], [1,4,3,2,5,6,7,8,9,10]});
        %   >>> H = G.setwiseStabilizer([2 3]);
        %   >>> H.order
        %       10
        %
        % Args:
        %   set (integer(1,\*)): The set to stabilize
        %
        % Returns:
        %   `+replab.PermutationGroup`: The subgroup that stabilizes the set
            s = replab.bsgs.SetwiseStabilizer(self, set).subgroup;
        end

        function s = pointwiseStabilizer(self, set)
        % Returns the subgroup that stabilizes the given set pointwise
        %
        % i.e. for this group ``G``, it returns ``H = {g \in G : g(i) = i, i \in set}``.
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> G = S4.pointwiseStabilizer([1 2]);
        %   >>> G.order == 2
        %       1
        %
        % Args:
        %   set (integer(1,\*)): The set to stabilize pointwise
        %
        % Returns:
        %   `+replab.PermutationGroup`: The subgroup that stabilizes the set pointwise
            c = self.chain.mutableCopy;
            c.baseChange(set, true);
            l = find(~ismember(c.base, set), 1); % find the first base point that is not in set
            if isempty(l)
                s = self.trivialSubgroup;
            else
                immutable = true;
                c = c.chainFromLevel(l, immutable);
                s = replab.PermutationGroup.fromChain(c, self.type);
            end
        end

        function o = orbits(self)
        % Returns the partition of the domain into orbits under this group
        %
        % Permutation group orbits are also called domains of transitivity,
        % see https://www.encyclopediaofmath.org/index.php/Transitive_group
        %
        % Returns:
        %   `.Partition`: The orbit partition
            G = zeros(self.nGenerators, self.domainSize);
            for i = 1:self.nGenerators
                G(i, :) = self.generators{i};
            end
            o = replab.Partition.permutationsOrbits(G);
        end

        function sub = stabilizer(self, p)
        % Returns the subgroup that stabilizes a given point
            sub = replab.PermutationGroup.fromChain(self.chain.stabilizer(p), self.type);
        end

        function w = wreathProduct(self, A)
        % Returns the wreath product of a compact group by this permutation group
        %
        % See https://en.wikipedia.org/wiki/Wreath_product
        %
        % Note that our notation is reversed compared to the Wikipedia page,
        % the permutation group is on the left hand side, as our convention
        % for semidirect product places the group acted upon on the right.
        %
        % Note that the return type depends on the argument type:
        % if ``A`` is a `.FiniteGroup`, the result will be a finite group too,
        % and if ``A`` is a `.NiceFiniteGroup`, the result will be of that type.
        %
        % Args:
        %   A (`.CompactGroup`): The group whose copies are acted upon
        %
        % Returns:
        %   `+replab.WreathProductGroup`: A wreath product group
            w = replab.WreathProductGroup.make(self, A);
        end

    end

    methods % Actions

        function A = naturalAction(self)
        % Returns the natural action of elements of this group on its domain
        %
        % This group natural domain is the set of integers ``{1..domainSize}``
        %
        % Returns:
        %   replab.Action: The natural action
            A = replab.perm.PermutationNaturalAction(self);
        end

        function A = vectorAction(self)
        % Returns the action of permutations on column vectors
        %
        % Acts on vectors of size `domainSize` by permuting their coefficients
        %
        % Returns:
        %   replab.Action: The vector action
            A = replab.perm.PermutationVectorAction(self);
        end

        function A = matrixAction(self)
        % Returns the simultaneous action of permutations on both rows and columns of square matrices
        %
        % Acts on matrices of size ``self.domainSize x self.domainSize``
        %
        % Returns:
        %   replab.Action: The matrix action
            A = replab.perm.PermutationMatrixAction(self);
        end

        function perm = indexRelabelingPermutation(self, g, indexRange)
        % Returns the permutation that acts by permuting tensor coefficients
        %
        % Let ``I = (i_1, ..., i_d)`` be a sequence of indices, where ``d = self.domainSize``
        % and ``1 <= i_1,...,i_d <= indexRange``
        %
        % We enumerate elements of ``I`` by first incrementing ``i_d``, then ``i_{d-1}``, etc...
        %
        % We compute the permutation of domain size ``indexRange^domainSize`` that acts on the
        % indices of I according to the argument ``g``.
        %
        % Args:
        %   g (permutation): Permutation of subindices
        %   indexRange (integer): Dimension of each subindex
        %
        % Returns:
        %   permutation: The permutation on the enumeration of indices
            n = self.domainSize;
            dims = indexRange * ones(1, n);
            perm = permute(reshape(1:prod(dims), dims), fliplr(n+1-g));
            perm = perm(:)';
        end

        function phi = indexRelabelingMorphism(self, indexRange)
        % Returns the morphism the permutation action of this group on tensor coefficients
        %
        % The tensor coefficients correspond to ``R^ir x R^ir ... (domainSize times)``
        % where ``ir = indexRange``.
        %
        % The enumeration of indices is done in the same order as in `.indexRelabelingPermutation`:
        % if ``I = (i1, ..., id)`` is a sequence of indices, we increment first ``id``, then
        % ``i_{d-1}`` and so on.
        %
        % See also:
        %   `+replab.PermutationGroup.indexRelabelingPermutation`
        %
        % Args:
        %   indexRange (integer): Dimension of each subindex
        %
        % Returns:
        %   function_handle: The permutation group homomorphism
            phi = @(g) self.indexRelabelingPermutation(g, indexRange);
        end

    end

    methods % Representations

        function rho = indexRelabelingRep(self, indexRange)
        % Representation that permutes the indices of a tensor
        %
        % It acts on the tensor space R^ir x R^ir ... (domainSize times)
        % where ir = indexRange, by permuting the indices.
        %
        % The representation returned is real.
        %
        % See also:
        %   `+replab.PermutationGroup.indexRelabelingPermutation`
        %
        % Args:
        %   indexRange (integer): Dimension of the tensor components/range of the subindices
        %
        % Returns:
        %   replab.Rep: The desired permutation representation
            rho = replab.rep.IndexRelabelingRep(self, indexRange);
        end

        function rho = naturalRep(self)
        % Returns the natural permutation representation of this permutation group
        %
        % Returns:
        %   replab.Rep: The (real) natural permutation representation
            rho = self.permutationRep(self.domainSize, 'preimages', self.generators, 'images', self.generators);
        end

        function rho = standardRep(self)
        % Returns the standard representation of this permutation group
        %
        % It is an abuse of terminology as the "standard representation" is
        % the faithful $n-1$ dimensional representation of the symmetric group
        % acting on $n$ elements; but we can reuse that subrepresentation on
        % subgroups of the symmetric group.
        %
        % It corresponds to the representation complementary to the
        % trivial representation with basis ``[1, 1, ..., 1]'``
        %
        % Returns:
        %   `+replab.SubRep`: The (real) standard representation
            d = self.domainSize;
            [injectiond, projection] = replab.sym.sageSpechtStandardBasis(d);
            if replab.init.cyclolab().works
                rho = self.naturalRep.subRep(replab.cyclotomic(injectiond)/d, 'projection', replab.cyclotomic(projection));
            else
                rho = self.naturalRep.subRep(injectiond/d, 'projection', projection);
            end
        end

        function rho = signRep(self)
        % Returns the sign representation of this permutation group
        %
        % Returns:
        %   `+replab.Rep`: One dimensional sign representation of this group
            preimages = self.generators;
            images = cellfun(@(g) replab.Permutation.sign(g), preimages, 'uniform', 0);
            rho = self.signedPermutationRep(1, 'preimages', preimages, 'images', images);
        end

    end

    methods (Static) % Construction of standard groups

        function G = symmetric(n)
        % Returns the symmetric group acting on a certain domain size
        %
        % Args:
        %   domainSize (integer): Domain size, must be >= 0
        %
        % Returns:
        %   `+replab.PermutationGroup`: Symmetric group
            G = replab.SymmetricGroup.make(n);
        end

        function G = kleinFourGroup
        % Constructs the Klein Four-Group
        %
        % This corresponds to the symmetry group of a non-square rectangle, and corresponds to the direct product ``S2 x S2``.
        %
        % Returns:
        %   `+replab.PermutationGroup`: The Klein four-group as a permutation gorup
            g1 = [2,1,4,3];
            g2 = [3,4,1,2];
            G = replab.PermutationGroup.of(g1, g2);
        end

        function G = dihedral(n)
        % Constructs the dihedral group of order ``2*n``
        %
        % This corresponds to the group of symmetries of the polygon with ``n`` vertices
        %
        % Args:
        %   n (integer): Half the dihedral group order
        %
        % Returns:
        %   `+replab.PermutationGroup`: The dihedral group permuting the vertices of the ``n``-gon
            switch n
              case 1
                G = replab.PermutationGroup.symmetric(2);
              case 2
                G = replab.PermutationGroup.kleinFourGroup;
              otherwise
                g1 = fliplr(1:n);
                g2 = [2:n 1];
                G = replab.PermutationGroup.of(g1, g2);
            end
        end

        function G = alternating(n)
        % Constructs the alternating group
        %
        % Args:
        %   n (integer): Group degree
        %
        % Returns:
        %   `+replab.PermutationGroup`: The alternating group of degree ``n``
            Sn = replab.PermutationGroup.symmetric(n);
            if n <= 2 % special case: group is trivial
                G = replab.PermutationGroup.trivial(n);
            else
                t = [2 3 1 4:n];
                if n == 3 % special case: it is a cyclic group, one generator only
                    G = Sn.subgroupWithGenerators({t}, Sn.order/2);
                    return
                end
                % generators from page 2100 of
                % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
                if mod(n, 2) == 0
                    s = [2 1 4:n 3];
                else
                    s = [1 2 4:n 3];
                end
                G = Sn.subgroupWithGenerators({s t}, Sn.order/2);
            end
        end

        function G = cyclic(n)
        % Constructs the cyclic group of order ``n`` acting on ``n`` points
        %
        % Args:
        %   n (integer): Cyclic group order and domain size
        %
        % Returns:
        %   `+replab.PermutationGroup`: The cyclic group of given order/domain size
            Sn = replab.PermutationGroup.symmetric(n);
            if n == 1
                G = Sn;
            else
                G = Sn.subgroupWithGenerators({[2:n 1]}, vpi(n));
            end
        end

        function G = trivial(n)
        % Constructs the trivial permutation group acting on ``n`` points
        %
        % Example:
        %   >>> G = replab.PermutationGroup.trivial(4);
        %   >>> G.order
        %     1
        %
        % Args:
        %   n (integer): Domain size
        %
        % Returns:
        %   `+replab.PermutationGroup`: Trivial group
            Sn = replab.S(n);
            G = Sn.subgroup({});
        end

    end

    methods (Static) % Construction of permutation groups

        function G = permutingGivenPoints(n, points)
        % Constructs the group that permutes the given points
        %
        % Essentially constructs the symmetric group of order ``|points|``.
        %
        % Args:
        %   n (integer): Domain size of the created group
        %   points (integer(1,\*)): Set of points being permuted
        %
        % Returns:
        %   `.PermutationGroup`: The permutation group
            switch length(points)
              case {0, 1}
                G = replab.PermutationGroup.trivial(n);
              case 2
                gen = 1:n;
                gen(points) = fliplr(points);
                G = replab.PermutationGroup.of(gen);
              otherwise
                gen1 = 1:n;
                gen2 = 1:n;
                gen1(points(1:2)) = gen1([points(2) points(1)]);
                gen2(points) = gen2([points(2:end) points(1)]);
                G = replab.PermutationGroup.of(gen1, gen2);
            end
        end

        function G = of(varargin)
        % Constructs a nontrivial permutation group from the given generators
        %
        % If you do not know the number of generators in advance, and would like to handle the
        % case of a trivial group, use ``Sn = replab.S(n); Sn.subgroup(generators)`` instead.
        %
        % Example:
        %   >>> G = replab.PermutationGroup.of([2 3 4 1], [4 3 2 1]);
        %   >>> G.order
        %     8
        %
        % Args:
        %   varargin (cell(1,\*) of permutation): Group generators
        %
        % Returns:
        %   `+replab.PermutationGroup`: The permutation group given as the closure of the generators
            assert(nargin > 0, 'Must be called with at least one generator');
            n = length(varargin{1});
            Sn = replab.S(n);
            G = Sn.subgroup(varargin);
        end

    end

end
