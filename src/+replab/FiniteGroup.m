classdef FiniteGroup < replab.CompactGroup & replab.FiniteSet
% Describes a group with a finite number of elements
%
% Each finite group has a type, that describes the most general group embedding its elements.
% For example, permutations of domain size ``n`` are embedded in the symmetric group of degree ``n``. We call that
% embedding a `.type`. All groups with the same type have the same `.niceMorphism` enabling the shift of computations
% to a "nice" permutation group.

    properties (SetAccess = protected)
        generators % (cell(1,\*) of `.type` elements): Group generators
    end

    methods % Implementations

        function res = ne(self, rhs)
            res = ~(self == rhs);
        end

        function res = eq(self, rhs)
            res = isa(rhs, 'replab.FiniteGroup') && self.hasSameTypeAs(rhs) && self.isSubgroupOf(rhs) && rhs.isSubgroupOf(self);
        end

        function res = isequal(self, rhs)
            res = self == rhs;
        end

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Group(self);
            names{1, end+1} = 'generators';
            names{1, end+1} = 'representative';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nGenerators
                names{1, end+1} = sprintf('generator(%d)', i);
                values{1, end+1} = self.generator(i);
            end
            r = self.fastRecognize;
            if ~isempty(r)
                names{1,end+1} = 'recognize';
                values{1,end+1} = r;
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.FiniteGroupLaws(self);
        end

        % Group

        function m = morphismByFunction(self, target, imageElementFun)
            imgs = cellfun(imageElementFun, self.generators, 'uniform', 0);
            m = self.morphismByImages(target, 'preimages', self.generators, 'images', imgs);
        end

        function m = isomorphismByFunctions(self, target, preimageElementFun, imageElementFun)
        % Constructs a group isomorphism using preimage/image functions
        %
        % Args:
        %   target (`replab.Group`): Target group
        %   preimageElementFun (function_handle): Returns the source element for a target element
        %   imageElementFun (function_handle): Returns the target element for a source element
            m = self.isomorphismByFunction(self, target, imageElementFun);
        end

        function m = innerAutomorphism(self, by)
        % Returns the morphism that corresponds to left conjugation by an element
        %
        % Args:
        %   by (element of `.type`): Element to conjugate the group with
        %
        % Returns:
        %   `+replab.FiniteMorphism`: Conjugation morphism
            generatorImages = cellfun(@(g) self.leftConjugate(by, g), self.generators, 'uniform', 0);
            m = self.isomorphismByImages(self, 'preimages', self.generators, 'images', generatorImages);
        end

    end

    methods (Access = protected)

        function o = computeOrder(self)
        % See `.order`
            error('Abstract');
        end

        function D = computeDecomposition(self)
        % See `.decomposition`
            error('Abstract');
        end

        function a = computeAbelianInvariants(self)
        % See `.abelianInvariants`
        %
        % Inspired by GAP System "AbelianInvariants"
            if self.isTrivial
                a = zeros(1, 0);
                return
            end
            a = zeros(1, 0);
            G = self;
            dsg = self.derivedSubgroup;
            primeDivisors = unique(double(factor(self.order)));
            for i = 1:length(primeDivisors)
                p = primeDivisors(i);
                ranks = zeros(1, 0);
                r = 0;
                while r ~= 1
                    H = dsg;
                    for j = 1:G.nGenerators
                        H = H.closure(G.composeN(G.generator(j), p));
                    end
                    r = G.order / H.order;
                    G = H;
                    if r ~= 1
                        ranks = [ranks length(factor(r))];
                    end
                end
                if ~isempty(ranks)
                    l = ones(1, ranks(1));
                    for i = ranks
                        l(1:i) = l(1:i) * p;
                    end
                    a = [a l];
                end
            end
            a = sort(a);
        end

        function e = computeExponent(self)
        % See `.exponent`
            eo = self.conjugacyClasses.classElementOrders;
            eo = unique(eo);
            e = eo(1);
            for i = 2:length(eo)
                e = lcm(e, eo(i));
            end
        end

        function c = computeCharacterTable(self)
        % See `.characterTable`
            r = self.recognize;
            assert(~isempty(r), 'This group has not been recognized.');
            c = r.atlasEntry.characterTable;
            assert(~isempty(c), 'This group does not have a stored character table.');
            c = c.imap(r.isomorphism);
        end

        function c = computeConjugacyClasses(self)
        % See `.conjugacyClasses`
            error('Abstract');
        end

        function R = computeFastRecognize(self)
            R = [];
            if self.niceMorphism.image.domainSize < replab.globals.fastChainDomainSize
                c = self.niceMorphism.image.partialChain;
                if ~c.isMutable
                    R = replab.Atlas.recognize(self);
                end
            end
        end

        function R = computeRecognize(self)
            R = replab.Atlas.recognize(self);
        end

        function res = computeIsCyclic(self)
            error('Abstract');
        end

        function res = computeIsSimple(self)
            error('Abstract');
        end

        function sub = computeDerivedSubgroup(self)
        % See `.derivedSubgroup`
            error('Abstract');
        end

    end

    methods % Group properties

        function s = nElements(self)
            s = self.order;
        end

        function o = order(self)
        % Returns the group order
        %
        % Returns:
        %   vpi: The group order
            o = self.cached('order', @() self.computeOrder);
        end

        function D = decomposition(self)
        % Returns a decomposition of this group as a product of sets
        %
        % Returns:
        %   `.FiniteGroupDecomposition`: The group decomposition
            D = self.cached('decomposition', @() self.computeDecomposition);
        end

        function e = exponent(self)
        % Returns the group exponent
        %
        % The group exponent is the smallest integer ``e`` such that ``g^e == identity`` for all ``g`` in ``G``.
        %
        % Returns:
        %   vpi: The group exponent
            e = self.cached('exponent', @() self.computeExponent);
        end

        function a = abelianInvariants(self)
        % Returns the group abelian invariants
        %
        % It computes the decomposition of the factor group of this group by its derived subgroup
        % (which is abelian), and returns the primary decomposition of that abelian group.
        %
        % Example:
        %   >>> G = replab.PermutationGroup.cyclic(100);
        %   >>> isequal(G.abelianInvariants, [4 25])
        %       1
        %
        % Returns:
        %   integer(1,\*): Group abelian invariants
            a = self.cached('abelianInvariants', @() self.computeAbelianInvariants);
        end

        function c = conjugacyClass(self, g)
        % Returns the conjugacy class corresponding to the given element
        %
        % Args:
        %   g (element): Arbitrary group element
        %
        % Returns:
        %   `.ConjugacyClass`: The conjugacy class containing the given element
            c = replab.ConjugacyClass.make(self, g);
        end

        function c = conjugacyClasses(self)
        % Returns the conjugacy classes of this group
        %
        % Returns:
        %   `+replab.ConjugacyClasses`: Conjugacy classes
            c = self.cached('conjugacyClasses', @() self.computeConjugacyClasses);
        end

        function c = characterTable(self)
            c = self.cached('characterTable', @() self.computeCharacterTable);
        end

        function R = fastRecognize(self)
        % Attempts to recognize this group in the standard atlas
        %
        % Returns:
        %   `+replab.AtlasResult` or []: A result in case the group is identified; or ``[]`` if unrecognized.
            R = self.cached('fastRecognize', @() self.computeFastRecognize);
        end

        function R = recognize(self)
        % Attempts to recognize this group in the standard atlas
        %
        % Returns:
        %   `+replab.AtlasResult` or []: A result in case the group is identified; or ``[]`` if unrecognized.
            R = self.cached('recognize', @() self.computeRecognize);
        end

        function b = isTrivial(self)
        % Tests whether this group is trivial
        %
        % Returns:
        %   logical: True if this group is trivial (i.e. has only one element)
            b = self.nGenerators == 0;
        end

        function res = isCommutative(self)
        % Returns whether this group is commutative
            for i = 1:self.nGenerators
                gi = self.generator(i);
                for j = 1:i-1
                    gj = self.generator(j);
                    if ~self.eqv(self.compose(gi, gj), self.compose(gj, gi))
                        res = false;
                        return
                    end
                end
            end
            res = true;
            return
        end

        function res = isCyclic(self)
        % Returns whether this group is a cyclic group
            res = self.cached('isCyclic', @() self.computeIsCyclic);
        end

        function res = isSimple(self)
        % Returns whether this group is simple
            res = self.cached('isSimple', @() self.computeIsSimple);
        end

    end

    methods % Group elements

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

        function l = factorizeLetters(self, element)
        % Factorizes an element as a word in the generators (returns letters)
        %
        % Args:
        %   element (element this group): Element to factorize
        %
        % Returns:
        %   integer(1,\*): Letters of the word in the generators
            error('Abstract');
        end

        function w = factorizeWord(self, element, names)
        % Factorizes an element as a word in the generators (returns word)
        %
        % Example:
        %   >>> G = replab.S(3);
        %   >>> G.factorizeWord([2 3 1])
        %       'x1'
        %   >>> names = {'s', 't'};
        %   >>> G.factorizeWord([2 3 1], names)
        %       's'
        %
        % Args:
        %   element (element of this group): Element to factorize
        %   names (cell(1,\*) of charstring, optional): Generator names, default `.defaultGeneratorNames`
        %
        % Returns:
        %   charstring: Word in the generators
            if nargin < 3 || isempty(names)
                names = self.defaultGeneratorNames;
            end
            l = self.factorizeLetters(element);
            w = replab.fp.Letters.print(l, names);
        end

        function g = imageWord(self, word, names)
        % Computes the image of a word in the group generators
        %
        % Example:
        %   >>> G = replab.S(3);
        %   >>> G.imageWord('x1')
        %       2     3     1
        %   >>> names = {'s', 't'};
        %   >>> G.imageWord('s', names)
        %       2     3     1
        %
        % Args:
        %   word (charstring): Word written using mathematical notation
        %   names (cell(1,\*) of charstring, optional): Generator names, default `.defaultGeneratorNames`
        %
        % Returns:
        %   element: Element of the group corresponding to the given word
            if nargin < 3 || isempty(names)
                names = self.defaultGeneratorNames;
            end
            l = replab.fp.Letters.parse(word, names);
            g = self.imageLetters(l);
        end

        function g = imageLetters(self, letters)
        % Returns the image of a word in the group generators, where the word is represented by letters
        %
        % The letters take the values ``1...nG`` which represent the group generators, and also the values
        % ``-nG...-1`` which represent the inverses of those generators.
        %
        % Args:
        %   letters (integer(1,\*)): Letters composing the word.
        %
        % Returns;
        %   element: Element of the group corresponding to the given word
            g = self.composeLetters(self.generators, letters);
        end

        function b = contains(self, g)
        % Tests whether this group contains the given element
        %
        % Args:
        %   g (element of `.type`): Element to test membership of
        %
        % Returns:
        %   logical: True if this group contains ``g`` and false otherwise
            error('Abstract');
        end

        function o = elementOrder(self, g)
        % Returns the order of a group element
        %
        % Args:
        %   g (element): Group element
        %
        % Returns:
        %   vpi: The order of ``g``, i.e. the smallest ``o`` such that ``g^o == identity``
            error('Abstract');
        end

    end

    methods % Construction of groups

        function conj = leftConjugateGroup(self, by)
        % Returns the left conjugate of the current group by the given element
        %
        % ``res = self.leftConjugateGroup(by)``
        %
        % In particular, it ensures that
        % ``res.generator(i) = self.type.leftConjugate(by, self.generator(i))``
        %
        % Args:
        %   by (element of `.type`): Element to conjugate the group with
        %
        % Returns:
        %   `+replab.FiniteGroup`: The conjugated group
            newGenerators = cellfun(@(g) self.type.leftConjugate(by, g), self.generators, 'uniform', 0);
            conj = self.type.subgroupWithGenerators(newGenerators, self.order);
        end

        function res = closure(self, obj)
        % Computes the group generated by the elements of this group and another object
        %
        % Example:
        %    >>> G = replab.PermutationGroup.cyclic(3);
        %    >>> S3 = replab.S(3);
        %    >>> H = S3.subgroup({[2 1 3]});
        %    >>> C = G.closure(H);
        %    >>> C == S3
        %        1
        %
        % Args:
        %   rhs (`.FiniteGroup` or group element): Element(s) to add
        %
        % Returns:
        %   `.FiniteGroup`: The closure
            error('Abstract');
        end

        function res = normalClosure(self, obj)
        % Computes the normal closure of an object in the closure of this group and this object
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> S4 = S5.subgroup({[2 3 4 1 5] [2 1 3 4 5]});
        %   >>> nc = S4.normalClosure(S5.subgroup({[1 2 4 5 3]}));
        %   >>> nc == S5.subgroup({[1 2 4 5 3] [5 2 3 1 4] [2 5 3 4 1]})
        %       1
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> nc = S4.normalClosure(S4.subgroup({[2 3 1 4]}));
        %   >>> nc == S4.subgroup({[2 3 1 4] [1 3 4 2]});
        %       1
        %
        % Args:
        %   obj (`+replab.FiniteGroup` or group element): Element to compute the normal closure with
        %
        % Returns:
        %   `+replab.FiniteGroup`: The normal closure of ``obj`` in this group
            error('Abstract');
        end

    end

    methods % Subgroups

        function sub = subgroupWithGenerators(self, generators, order)
        % Constructs a subgroup of the current group from generators
        %
        % Guarantees that the constructed subgroup preserves the order of the generators. However, the
        % given generators cannot contain the identity element.
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> G = S5.subgroupWithGenerators({[2 3 4 5 1]});
        %   >>> isequal(G.generators{1}, [2 3 4 5 1])
        %       1
        %
        % Args:
        %   generators (cell(1,\*) of elements of this group): List of generators, does not contain the identity
        %   order (vpi or ``[]``, optional): Subgroup order
        %
        % Returns:
        %   `+replab.FiniteGroup`: The subgroup generated by `.generators`
            error('Abstract');
        end

        function sub = subgroup(self, elements, order)
        % Constructs a subgroup of the current group from elements
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> G = S5.subgroup({[1 2 3 4 5]});
        %   >>> isempty(G.generators)
        %       1
        %
        % Args:
        %   elements (cell(1,\*) of elements of this group): List of elements
        %   order (vpi, optional): Subgroup order
        %
        % Returns:
        %   `+replab.FiniteGroup`: The subgroup generated by `.generators`
            if nargin < 3
                order = [];
            end
            mask = cellfun(@(e) self.isIdentity(e), elements);
            sub = self.subgroupWithGenerators(elements(~mask), order);
        end

        function sub = randomProperSubgroup(self, nSteps)
        % Constructs a random proper subgroup of this group
        %
        % Args:
        %   nSteps (integer, optional): How many steps of reduction to perform, default 1
        %
        % Returns:
        %   `.FiniteGroup`
            if nargin < 2 || isempty(nSteps)
                nSteps = 1;
            end
            sub = self;
            for i = 1:nSteps
                if isprime(sub.order)
                    if i == 1
                        error('This cyclic group of prime order has no proper subgroup');
                    else
                        return
                    end
                end
                candidate = sub;
                while candidate.order == sub.order
                    s1 = sub.sample;
                    while sub.isIdentity(s1)
                        s1 = sub.sample;
                    end
                    s2 = sub.sample;
                    while sub.isIdentity(s2)
                        s2 = sub.sample;
                    end
                    candidate = sub.subgroup({s1 s2});
                end
                sub = candidate;
            end
        end

        function sub = trivialSubgroup(self)
        % Returns the trivial subgroup of this group
        %
        % Returns:
        %   `+replab.FiniteGroup`: The trivial subgroup
            sub = self.subgroup({}, vpi(1));
        end

        function sub = derivedSubgroup(self)
        % Computes the derived subgroup of this group
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> S4.order
        %       24
        %   >>> D = S4.derivedSubgroup;
        %   >>> D == replab.PermutationGroup.alternating(4)
        %       1
        % Returns:
        %   `+replab.FiniteGroup`: The derived subgroup
            sub = self.cached('derivedSubgroup', @() self.computeDerivedSubgroup);
        end

        function ser = derivedSeries(self)
            D = self.derivedSubgroup;
            if self == D
                ser = {self};
            else
                ser = horzcat({self}, D.derivedSeries);
            end
        end

        function sub = center(self)
        % Returns the center of this group
        %
        % Returns:
        %   `.FiniteGroup`: The center
            sub = self.cached('center', @() self.centralizer(self));
        end

        function sub = centralizer(self, arg)
        % Returns the centralizer of the given object in this group
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizer([2 3 1 4]);
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizer(G.subgroup({[2 3 1 4]}));
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Args:
        %   rhs (`.FiniteGroup` or group element): Element to compute the centralizer of
        %
        % Returns:
        %   `.FiniteGroup`: The centralizer
            error('Abstract');
        end

        function sub = intersection(self, G)
        % Computes the intersection of two groups
        %
        % Example:
        %   >>> D12 = replab.S(6).subgroup({[1 6 5 4 3 2], [2 1 6 5 4 3]});
        %   >>> G = replab.S(6).subgroup({[2 1 3 4 5 6], [2 3 4 5 1 6]});
        %   >>> I = D12.intersection(G);
        %   >>> I == replab.S(6).subgroup({[5 4 3 2 1 6]})
        %       1
        %
        % Args:
        %   G (`+replab.FiniteGroup`): Group with the same type as this group
        %
        % Returns:
        %   `+replab.FiniteGroup`: Subgroup representing the intersection
            error('Abstract');
        end

    end

    methods % Cosets

        function l = isNormalizedBy(self, element)
        % Returns whether a given element/group normalizes this group
        %
        % This is true when ``element * group * element^-1 == group``.
        %
        % The same definition when ``element`` is a group; then we ask that all this group elements satisfy that property.
        %
        % Args:
        %   element (group element or `.FiniteGroup`): Group element or group
        %
        % Returns:
        %   logical: True if the given element/group normalizes this group
            if isa(element, 'replab.FiniteGroup')
                l = all(cellfun(@(g) self.isNormalizedBy(g), element.generators));
            else
                l = false;
                for i = 1:self.nGenerators
                    g = self.generator(i);
                    if ~self.contains(self.type.leftConjugate(element, g))
                        return
                    end
                end
                l = true;
            end
        end

        function c = doubleCosets(self, H, K)
        % Returns the set of double cosets in this group by the given groups
        %
        % Args:
        %   H (`+replab.FiniteGroup`): First subgroup
        %   K (`+replab.FiniteGroup`): Second subgroup
        %
        % Returns:
        %   `+replab.DoubleCosets`: The set of double cosets
            c = replab.DoubleCosets(self, H, K);
        end

        function c = doubleCoset(self, element, K, parent)
        % Returns a double coset
        %
        % Returns the set ``self * element * K``.
        %
        % Args:
        %   element (group element): Group element
        %   K (`.FiniteGroup`): Finite group
        %   parent (`.FiniteGroup`, optional): Group containing all of ``self``, ``element`` and ``K``
        %
        % Returns:
        %   `+replab.DoubleCoset`: The constructed double coset
            if nargin < 4 || isempty(parent)
                parent = [];
            end
            c = replab.DoubleCoset.make(self, element, K, parent);
        end

        function c = normalCoset(self, element, parent)
        % Returns a normal coset
        %
        % Returns the set ``element * self == self * element``.
        %
        % Args:
        %   element (group element): Group element
        %   parent (`.FiniteGroup`, optional): Group containing both ``self`` and ``element``
        %
        % Returns:
        %   `+replab.NormalCoset`: The constructed normal coset
            assert(self.isNormalizedBy(element), 'The given element does not define a normal coset');
            if nargin < 3 || isempty(parent)
                parent = [];
            end
            c = replab.NormalCoset.make(self, element, parent);
        end

        function c = normalCosetsOf(self, subgroup)
        % Returns the set of normal cosets of the given subgroup in this group
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Normal subgroup of this group
        %
        % Returns:
        %   `.NormalCosets`: The set of normal cosets
            assert(subgroup.isNormalSubgroupOf(self), 'The given subgroup must be normal in parent group');
            c = replab.NormalCosets(self, subgroup);
        end

        function c = rightCoset(self, element, parent)
        % Returns a right coset
        %
        % Returns the set ``self * element``.
        %
        % Args:
        %   element (group element): Group element
        %   parent (`.FiniteGroup`, optional): Group containing both ``self`` and ``element``
        %
        % Returns:
        %   `+replab.RightCoset`: The constructed right coset
            if nargin < 3 || isempty(parent)
                parent = [];
            end
            c = replab.RightCoset.make(self, element, parent);
        end

        function c = rightCosetsOf(self, subgroup)
        % Returns the set of right cosets of the given subgroup in this group
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Subgroup of this group
        %
        % Returns:
        %   `+replab.RightCosets`: Right cosets
            c = replab.RightCosets(self, subgroup);
        end

        function c = mldivide(self, supergroup)
        % Shorthand for `.rightCosetsOf`
            c = supergroup.rightCosetsOf(self);
        end

        function c = leftCoset(self, element, parent)
        % Returns a left coset
        %
        % Returns the set ``element * self``.
        %
        % Args:
        %   element (group element): Group element
        %   parent (`.FiniteGroup`, optional): Group containing both ``self`` and ``element``
        %
        % Returns:
        %   `+replab.LeftCoset`: The constructed right coset
            if nargin < 3 || isempty(parent)
                parent = [];
            end
            c = replab.LeftCoset.make(self, element, parent);
        end

        function c = leftCosetsOf(self, subgroup)
        % Returns the set of left cosets of the given subgroup in this group
        %
        % Args:
        %   subgroup (`+replab.FiniteGroup`): Subgroup of this group
        %
        % Returns:
        %   `+replab.LeftCosets`: Left cosets
            c = replab.LeftCosets(self, subgroup);
        end

        function c = mrdivide(self, subgroup)
        % Shorthand for `.leftCosetsOf`
            c = self.leftCosetsOf(subgroup);
        end

        function B = findLeftConjugations(self, s, t, sCentralizer, tCentralizer)
        % Returns the set of all elements that left conjugates an element to another element
        %
        % Let ``s`` and ``t`` be two elements of this group. This returns the set of all elements
        % ``b`` such that ``t = b s b^-1`` or ``t = leftConjugate(b, s)``.
        %
        % When no such ``b`` exists, this returns ``[]``.
        %
        % Args:
        %   s (group element): Source element
        %   t (group element): Target element
        %   sCentralizer (`+replab.FiniteGroup` or ``[]``, optional): Centralizer of ``s`` in this group
        %   tCentralizer (`+replab.FiniteGroup` or ``[]``, optional): Centralizer of ``t`` in this group
        %
        % Returns:
        %   `+replab.LeftCoset` or ``[]``: Set of all elements of this group that left conjugates ``s`` to ``t`` if it exists
            error('Abstract');
        end

    end

    methods % Relations to other groups

        function res = isSubgroupOf(self, rhs)
        % Returns whether this group is a subgroup of another group
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> G = S3.subgroup({[2 1 3]});
        %   >>> G.isSubgroupOf(S3)
        %       1
        %
        % Args:
        %   rhs (`+replab.FiniteGroup`): Other group with the same type as this one
        %
        % Returns:
        %   logical: True if this group is a subgroup of ``rhs``
            res = all(cellfun(@(g) rhs.contains(g), self.generators));
        end

        function res = isNormalSubgroupOf(self, rhs)
        % Returns whether this group is a normal subgroup of another group
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> A3 = replab.PermutationGroup.alternating(3);
        %   >>> A3.isNormalSubgroupOf(S3)
        %       1
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> G = S3.subgroup({[2 1 3]});
        %   >>> G.isNormalSubgroupOf(S3)
        %       0
        %
        % Args:
        %   rhs (`+replab.FiniteGroup`): Other group with the same type as this one
        %
        % Returns:
        %   logical: True if this group is a normal subgroup of ``rhs``
            if ~self.isSubgroupOf(rhs)
                res = false;
                return
            end
            for i = 1:self.nGenerators
                subi = self.generator(i);
                for j = 1:rhs.nGenerators
                    rhsj = rhs.generator(j);
                    if ~self.contains(self.type.leftConjugate(rhsj, subi))
                        res = false;
                        return
                    end
                end
            end
            res = true;
        end

    end

    methods (Access = protected)

        function G = computeNiceGroup(self)
        % See `.niceGroup`
            error('Abstract');
        end

        function m = computeNiceMorphism(self)
        % See `.niceMorphism`
            error('Abstract');
        end

        function A = computeDefaultAbstractGroup(self)
        % See `.abstractGroup`
            error('Abstract');
        end

        function m = computeDefaultAbstractMorphism(self)
        % See `.abstractMorphism`
            error('Abstract');
        end

    end

    methods % Isomorphic groups

        function names = defaultGeneratorNames(self)
        % Returns default generator names
        %
        % Returns:
        %   cell(1,\*) of charstring: A list of generator names starting with ``x1``, ``x2``, ...
            names = arrayfun(@(i) ['x' num2str(i)], 1:self.nGenerators, 'uniform', 0);
        end

        function G = niceGroup(self)
        % Returns the permutation group isomorphic to this finite group
        %
        % The generators of the nice group must be in one-to-one correspondance with the generators of this group.
        %
        % This is the image of `.niceMorphism`.
        %
        % Returns:
        %   `+replab.PermutationGroup`: A permutation group isomorphic to this group
            G = self.cached('niceGroup', @() self.computeNiceGroup);
        end

        function G = abstractGroup(self, names)
        % Returns an abstract group isomorphic to this finite gorup
        %
        % The generators of the abstract group must be in one-to-one correspondance with the generators of this group.
        %
        % This is the image of `.abstractMorphism`; and the `+replab.AbstractGroup.niceGroup` will be the same as this group
        % `.niceGroup`.
        %
        % Args:
        %   names (cell(1,\*) of charstring or ``[]``, optional): Generator names, default to ``x1`` ... ``xN``.
        %
        % Returns:
        %   `+replab.AbstractGroup`: An abstract group isomorphic to this group
            if nargin < 2 || isempty(names)
                G = self.cached('defaultAbstractGroup', @() self.computeDefaultAbstractGroup);
            else
                G = self.abstractGroup.withRenamedGenerators(names);
            end
        end

    end


    methods % Morphisms

        function res = findIsomorphism(self, to)
        % Finds an isomorphism from this finite group to another finite group, if it exists
        %
        % Args:
        %   to (`+replab.FiniteGroup`): Target of the isomorphism
        %
        % Returns:
        %   `.FiniteIsomorphism` or ``[]``: The isomorphism if it exists, or an empty array
            res = self.findIsomorphisms(to, 'single', true, 'upToConjugation', true);
            if ~isempty(res)
                res = res{1};
            else
                res = [];
            end
        end

        function res = findIsomorphisms(self, to, varargin)
        % Finds all the isomorphisms from this finite group to another finite group
        %
        % Example:
        %   >>> G = replab.S(6);
        %   >>> m = G.findIsomorphisms(G, 'upToConjugation', true);
        %   >>> length(m)
        %       2
        %
        % Args:
        %   to (`+replab.FiniteGroup`): Target of the isomorphism
        %
        % Keyword Args:
        %   upToConjugation (logical, optional): Whether to list morphisms up to conjugation of the image group, default: false
        %   single (logical, optional): Whether to return maximum a single result, default: false
        %
        % Returns:
        %   cell(1,\*) of `.FiniteIsomorphism`: The isomorphisms
            args = struct('upToConjugation', false, 'single', false);
            args = replab.util.populateStruct(args, varargin);
            if args.single
                args.upToConjugation = true;
            end
            F = self.abstractGroup;
            G = to;
            A = to;
            if F.order ~= G.order
                res = cell(1, 0);
                return
            end
            fm = replab.mrp.FindMorphisms(F, G, A, 'isomorphisms', args.single);
            if args.upToConjugation || args.single
                res = fm.searchUpToConjugation;
            else
                res = fm.searchAll;
            end
            res = cellfun(@(m) self.abstractMorphism.andThen(m), res, 'uniform', 0);
            res = cellfun(@(m) m.toIsomorphism, res, 'uniform', 0);
        end

        function res = findMorphisms(self, to, varargin)
        % Finds all the morphisms from this finite group to another finite group
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> S4 = replab.S(4);
        %   >>> m = S4.findMorphisms(S3, 'upToConjugation', true, 'surjective', false);
        %   >>> length(m)
        %       3
        %
        % Args:
        %   to (`+replab.FiniteGroup`): Target of the morphism
        %
        % Keyword Args:
        %   upToConjugation (logical, optional): Whether to list morphisms up to conjugation of the image group, default: false
        %   surjective (logical, optional): Whether to consider only surjective morphisms (or epimorphisms), whose image span ``to``, default: false
        %   single (logical, optional): Whether to return maximum a single result, default: false
        %
        % Returns:
        %   cell(1,\*) of `.FiniteMorphism`: The morphisms
            args = struct('upToConjugation', false, 'surjective', false, 'single', false);
            args = replab.util.populateStruct(args, varargin);
            if args.single
                args.upToConjugation = true;
            end
            if args.surjective
                if self.order == to.order
                    res = self.findIsomorphisms(to, 'upToConjugation', args.upToConjugation, 'single', args.single);
                    return
                else
                    filter = 'epimorphisms';
                end
            else
                filter = 'morphisms';
            end
            F = self.abstractGroup;
            G = to;
            A = to;
            fm = replab.mrp.FindMorphisms(F, G, A, filter, args.single);
            if args.upToConjugation || args.single
                res = fm.searchUpToConjugation;
            else
                res = fm.searchAll;
            end
            res = cellfun(@(m) self.abstractMorphism.andThen(m), res, 'uniform', 0);
        end


        function m = conjugatingAutomorphism(self, by)
        % Returns the morphism that corresponds to left conjugation by an element
        %
        % Args:
        %   by (element of `.type`): Element to conjugate the group with
        %
        % Returns:
        %   `+replab.FiniteMorphism`: Conjugating automorphism
            generatorImages = cellfun(@(g) self.type.leftConjugate(by, g), self.generators, 'uniform', 0);
            assert(all(cellfun(@(g) self.contains(g), generatorImages)));
            m = self.morphismByImages(self, 'preimages', self.generators, 'images', generatorImages, 'nChecks', 0);
        end

        function f = niceMorphism(self)
        % Returns the isomorphism from this group to a permutation group
            f = self.cached('niceMorphism', @() self.computeNiceMorphism);
        end

        function m = abstractMorphism(self, names)
        % Returns an isomorphism to an abstract group
        %
        % Example:
        %   >>> G = replab.S(3);
        %   >>> f = G.abstractMorphism({'s' 't'});
        %   >>> f.imageElement([2 3 1])
        %       's'
        %
        % Args:
        %   names (cell(1,\*) of charstring, optional): Generator names, default `.defaultGeneratorNames`
        %
        % Returns:
        %   `.FiniteIsomorphism`: Isomorphism to an abstract group
            if nargin < 2 || isempty(names)
                m = self.cached('defaultAbstractGroupIsomorphism', @() self.computeDefaultAbstractMorphism);
            else
                A = self.abstractGroup(names);
                m = self.niceMorphism.andThen(A.niceMorphism.inverse);
            end
        end

        function l = isMorphismByImages(self, target, varargin)
        % Checks whether the given images describe a group morphism
        %
        % The calling convention of this method is the same as `.morphismByImages`, except that ``nChecks`` is missing.
        %
        % Args:
        %   target (`+replab.Group`): Target group
        %
        % Keyword Args:
        %   preimages (cell(1, \*) of ``self`` elements): Preimages of the morphism which generate ``self``, defaults to ``self.generators``
        %   images (cell(1, \*) of ``target`` elements): Images of the given preimages, defaults to ``target.generators`` if ``target`` is a `.FiniteGroup`
        %
        % Returns:
        %   logical: True if the given images define a morphism
            assert(isa(target, 'replab.Group'));
            args = struct('preimages', {self.generators});
            if isa(target, 'replab.FiniteGroup')
                args.images = target.generators;
            else
                args.images = [];
            end
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.images) && ~self.isTrivial
                error('Images must be provided');
            end
            assert(length(args.preimages) == length(args.images), 'Number of images does not match the number of preimages');
            preId = cellfun(@(g) self.isIdentity(g), args.preimages);
            l = cellfun(@(g) target.isIdentity(g), args.images(preId));
            if ~all(l)
                return
            end
            l = self.isMorphismByImages_(target, args.preimages(~preId), args.images(~preId));
        end

        function m = morphismByImages(self, target, varargin)
        % Constructs a morphism to a group using images of generators
        %
        % The type of the morphism depends on the type of ``target``:
        %
        % * if ``target`` is a `.FiniteGroup`, this method returns a `.FiniteMorphism`,
        %
        % * if ``target`` is only a `.Group`, this method returns a `.Morphism`.
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> m = S4.morphismByImages(replab.S(3), 'images', {[1 3 2] [3 2 1]});
        %   >>> m.laws.checkSilent
        %       1
        %
        % Args:
        %   target (`.Group` or `.FiniteGroup`): Target of the morphism, the morphism image is a subgroup of this
        %
        % Keyword Args:
        %   preimages (cell(1, \*) of ``self`` elements): Preimages of the morphism which generate ``self``, defaults to ``self.generators``
        %   images (cell(1, \*) of ``target`` elements): Images of the given preimages, defaults to ``target.generators`` if ``target`` is a `.FiniteGroup`
        %   nChecks (integer or ``inf``): Number of randomized image checks to perform, if ``inf`` computes and verifies a presentation of ``self``
        %
        % Returns:
        %   `.Morphism` or `.FiniteMorphism`: The constructed morphism
            args = struct('nChecks', replab.globals.morphismNChecks, 'preimages', {self.generators});
            if length(varargin) == 1 && iscell(varargin{1})
                warning('Old call style deprecated, add a ''images'' keyword argument');
                args.images = varargin{1};
            else
                if isa(target, 'replab.FiniteGroup')
                    args.images = target.generators;
                else
                    args.images = [];
                end
                args = replab.util.populateStruct(args, varargin);
            end
            if isempty(args.images) && ~self.isTrivial
                error('Images must be provided');
            end
            assert(length(args.preimages) == length(args.images), 'Number of images does not match the number of preimages');
            preId = cellfun(@(g) self.isIdentity(g), args.preimages);
            assert(all(cellfun(@(g) target.isIdentity(g), args.images(preId))), 'Images of identity must be identity');
            if isinf(args.nChecks)
                assert(self.isMorphismByImages_(target, args.preimages, args.images), 'The given images do not define a morphism');
            end
            m = self.morphismByImages_(target, args.preimages(~preId), args.images(~preId));
            if isfinite(args.nChecks) && args.nChecks > 0
                for i = 1:args.nChecks
                    s1 = self.sample;
                    s2 = self.sample;
                    s12 = self.compose(s1, s2);
                    t1 = m.imageElement(s1);
                    t2 = m.imageElement(s2);
                    t12 = m.imageElement(s12);
                    t1_2 = target.compose(t1, t2);
                    assert(target.eqv(t12, t1_2), 'The given images do not define a morphism');
                end
            end
        end

        function m = isomorphismByImages(self, target, varargin)
        % Constructs an isomorphism to a group using images of generators
            m = self.morphismByImages(target, varargin{:}).toIsomorphism;
        end

        function m = isomorphismByFunction(self, target, imageElementFun)
        % Constructs an isomorphism to a group using an image function
            imgs = cellfun(imageElementFun, self.generators, 'uniform', 0);
            m = self.isomorphismByImages(target, 'preimages', self.generators, 'images', imgs);
        end

        function g1 = imap(self, f)
        % Maps this finite group through an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.isSubgroupOf(f.source)``
            g1 = f.imageGroup(self);
        end

    end

    methods (Access = protected)

        function m = morphismByImages_(self, target, preimages, images)
        % Implements the `.morphismByImages` method
        %
        % Does not perform checks; preimages must not contain the identity
            error('Abstract');
        end

        function l = isMorphismByImages_(self, target, preimages, images)
        % Implements the `.isMorphismByImages` method
        %
        % Does not perform checks; preimages must not contain the identity
            source = self.subgroupWithGenerators(preimages);
            assert(source.order == self.order, 'The morphism preimages do not generate the source group');
            l = source.abstractGroup.isMorphismByImages(target, 'images', images);
        end

    end


    methods % Representations

        function rep = regularRep(self)
        % Returns the left regular representation of this group
        %
        % Returns:
        %   `.Rep`: The left regular representation as a real permutation representation
            error('Abstract');
        end

        function rho = repByImages(self, field, dimension, varargin)
        % Constructs a finite dimensional representation of this group from preimages/images pairs
        %
        % The preimages need to generate this group. If the preimages are omitted, they default to
        % the generators of this group.
        %
        % The images can either be exact, in which case they should be provided as `.cyclotomic` matrices,
        % or as ``double`` matrices -- the latter only if they have integer entries.
        % In the case of exact images, the keyword argument ``imagesErrorBound`` can be omitted.
        % When images are exact, RepLAB has the freedom of performing matrix products containing an
        % arbitrary number of factors without loss of precision, and thus a wider range of algorithm
        % can be used. Also, the constructed representation can return exact images using ``rep.image(g, 'exact')``
        %
        % Otherwise, the provided images will be considered inexact. In that case, RepLAB uses a factorization
        % algorithm to decompose group elements as short words in the preimages. The keyword argument
        % ``imagesErrorBound`` should be given as well; it corresponds to an upper bound on the error of the
        % matrix images in the Frobenius norm.
        % Either a single number can be provided, as an upper bound on all images in the provided set; or an individual
        % bound for each image. However, this error bound will be ignored for each image provided as an interval matrix.
        % If inexact floating-point images are provided but ``imagesErrorBound`` is omitted, an ad-hoc error is
        % estimated and a warning is emitted.
        %
        % If the images are exact, the ``isUnitary`` keyword parameter can be omitted.
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> m = S4.repByImages('R', 1, 'images', {-1 -1});
        %   >>> m.laws.checkSilent
        %       1
        %
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %
        % Keyword Args:
        %   preimages (cell(1,n) of ``self`` elements, optional): Preimages of the representation map which generate ``self``, defaults to ``self.generators``
        %   images (cell(1,n) of double(d,d) or cyclotomic(d,d), may be sparse): Images of the given preimages
        %   imagesErrorBound (double or double(1,d) or ``[]``): Error
        %   isUnitary (logical or ``[]``, optional): Value of the constructed `.Rep.isUnitary`; a "true" value accelerates computations
        %   isIrreducible (logical or ``[]``, optional): Value of the constructed `.Rep.isIrreducible`
        %   trivialDimension (integer or ``[]``, optional): Value of the constructed `.Rep.trivialDimension`
        %   frobeniusSchurIndicator (integer or ``[]``, optional): Value of the constructed `.Rep.frobeniusSchurIndicator`
        %   isDivisionAlgebraCanonical (logical or ``[]``, optional): Value of the constructed `.Rep.isDivisionAlgebraCanonical`
        %
        % Returns:
        %   `+replab.RepByImages`: The constructed group representation
            args = struct('preimages', {self.generators}, 'images', {{}});
            if length(varargin) == 1 && iscell(varargin{1})
                warning('Deprecated call convention');
                rho = replab.rep.repByImages(self, field, dimension, 'images', varargin{1});
            else
                rho = replab.rep.repByImages(self, field, dimension, varargin{:});
            end
        end

        function rho = permutationRep(self, dimension, varargin)
        % Constructs a permutation representation of this group
        %
        % The returned representation is real. Use `+replab.Rep.complexification` to obtain a complex representation.
        %
        % Args:
        %   dimension (integer): Dimension of the representation
        %
        % Keyword Args:
        %   preimages (cell(1,n) of ``self`` elements, optional): Preimages of the representation map which generate ``self``, defaults to ``self.generators``
        %   images (cell(1,\*) of permutations): Images of the generators as permutations of size ``dimension``
        %
        % Returns:
        %   `+replab.Rep`: The constructed group representation
            args = struct('preimages', {self.generators}, 'images', {{}});
            args = replab.util.populateStruct(args, varargin);
            preimages = args.preimages;
            images = cellfun(@(g) replab.Permutation.toSparseMatrix(g), args.images, 'uniform', 0);
            assert(length(preimages) == length(images), 'Must provide as many images as preimages');
            imageElements = cellfun(@(g) [g; zeros(1, dimension)], args.images, 'uniform', 0);
            imageGroup = replab.perm.GeneralizedSymmetricGroup(dimension, 1);
            rho = replab.rep.RepByImages_monomial(self, 'R', dimension, preimages, images, imageGroup, imageElements, ...
                                                  'isIrreducible', dimension == 1);
        end

        function rho = signedPermutationRep(self, dimension, varargin)
        % Returns a real signed permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification`` to obtain a complex representation.
        %
        % Args:
        %   dimension: Dimension of the representation
        %
        % Keyword Args:
        %   preimages (cell(1,n) of ``self`` elements, optional): Preimages of the representation map which generate ``self``, defaults to ``self.generators``
        %   images (cell(1,\*) of signed permutations): Images of the generators as signed permutations of size ``dimension``
        %
        % Returns:
        %   `+replab.Rep`: The constructed group representation
            args = struct('preimages', {self.generators}, 'images', {{}});
            args = replab.util.populateStruct(args, varargin);
            preimages = args.preimages;
            images = cellfun(@(g) replab.SignedPermutation.toSparseMatrix(g), args.images, 'uniform', 0);
            assert(length(preimages) == length(images), 'Must provide as many images as preimages');
            imageElements = cellfun(@(g) [abs(g); (1-sign(g))/2], args.images, 'uniform', 0);
            imageGroup = replab.perm.GeneralizedSymmetricGroup(dimension, 2);
            rho = replab.rep.RepByImages_monomial(self, 'R', dimension, preimages, images, imageGroup, imageElements);
        end

    end

end
