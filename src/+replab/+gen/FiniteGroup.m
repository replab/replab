classdef FiniteGroup < replab.FiniteGroup & replab.gen.FiniteSet
% A generic finite group

    methods

        function self = FiniteGroup(type, generators, varargin)
        % Constructs a nice finite group
        %
        % If ``nice`` is provided, then ``generatorNames``, ``order``, ``relators`` cannot be provided.
        % If ``nice`` is provided, then ``niceIsomorphism`` must be provided too.
        %
        % Args:
        %   type (`+replab.FiniteGroupType`): Finite group type
        %   generators (cell(1,\*) of group elements): Generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring, optional): Names of the generators
        %   nice (`+replab.FiniteGroup`): Nice group on which computations are performed
        %   niceIsomorphism (`+replab.+gen.NiceIsomorphism`): Isomorphism
        %   order (vpi or integer, optional): Group order
        %   relators (cell(1,\*) of charstring or integer(1,\*), optional): Relators
            args = struct('generatorNames', [], 'nice', [], 'niceIsomorphism', [], 'order', [], 'relators', []);
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.nice)
                if isempty(args.niceIsomorphism)
                    niceIsomorphism = type.niceIsomorphism(generators);
                else
                    niceIsomorphism = args.niceIsomorphism;
                end
                niceArgs = cell(1, 0);
                if ~isempty(args.generatorNames)
                    niceArgs{1,end+1} = 'generatorNames';
                    niceArgs{1,end+1} = args.generatorNames;
                end
                if ~isempty(args.order)
                    niceArgs{1,end+1} = 'order';
                    niceArgs{1,end+1} = args.order;
                end
                if ~isempty(args.relators)
                    niceArgs{1,end+1} = 'relators';
                    niceArgs{1,end+1} = args.relators;
                end
                niceGenerators = cellfun(@(g) niceIsomorphism.imageElement(g), generators, 'uniform', 0);
                nice = niceIsomorphism.target.type.groupWithGenerators(niceGenerators, niceArgs{:});
            else
                assert(~isempty(args.niceIsomorphism));
                assert(isempty(args.generatorNames));
                assert(isempty(args.order));
                assert(isempty(args.relators));
                nice = args.nice;
                niceIsomorphism = args.niceIsomorphism;
            end
            self.type = type;
            self.nice = nice;
            self.niceIsomorphism = niceIsomorphism;
            self.identity = type.identity;
            self.representative_ = type.identity;
            self.generators = generators;
        end

    end

    methods (Access = protected) % Implementations

        % Morphisms

        function m = morphismByImages_(self, target, preimages, images, imageElementFun)
            first = self.niceIsomorphism.restrictedSource(self); % maps this to the perm group
            preimagesNG = cellfun(@(g) first.imageElement(g), preimages, 'uniform', 0);
            second = self.nice.morphismByImages(target, 'preimages', preimagesNG, 'images', images, 'nChecks', 0); % from the perm group to the images
            if isa(second, 'replab.FiniteMorphism')
                m = replab.mrp.FiniteComposition(second, first, imageElementFun);
            else
                m = replab.mrp.Composition(second, first, imageElementFun);
            end
        end

    end

    methods % Implementation

        % Str

        function s = headerStr(self)
            if self.knownOrder
                s = sprintf('Finite group of order %s', strtrim(num2str(self.order)));
            else
                s = 'Finite group';
            end
        end

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteGroup(self);
        end

        function names = hiddenFields(self)
            names1 = hiddenFields@replab.gen.FiniteSet(self);
            names2 = hiddenFields@replab.FiniteGroup(self);
            names = union(names1, names2);
        end

        % Obj

        function l = laws(self)
            l = laws@replab.FiniteGroup(self);
        end

        % Domain

        function b = eqv(self, x, y)
            b = self.type.eqv(x, y);
        end

        function g = sample(self)
            g = self.niceIsomorphism.preimageElement(self.nice.sample);
        end

        % Monoid

        function z = compose(self, x, y)
            z = self.type.compose(x, y);
        end

        % Group

        function y = inverse(self, x)
            y = self.type.inverse(x);
        end

        % gen.FiniteSet

        function l = compatibleWithNiceIsomorphism(self, iso)
            l = false; % guilty until proven innocent
            if ~iso.target.hasSameTypeAs(self.nice)
                return
            end
            type = self.nice.type;
            for i = 1:self.nice.nGenerators
                g = self.generator(i);
                if ~iso.sourceContains(g)
                    return
                end
                imgGen = iso.imageElement(g);
                niceGen = self.nice.generator(i);
                if ~type.eqv(imgGen, niceGen)
                    return
                end
            end
            l = true;
        end

        % FiniteSet

        function g1 = imap(self, f)
            g1 = imap@replab.FiniteGroup(self, f);
        end

        % FiniteGroup

        function a = abelianInvariants(self)
            a = self.nice.abelianInvariants;
        end

        function c1 = centralizer(self, obj)
            if isa(obj, 'replab.FiniteGroup')
                c = self.nice.centralizer(self.niceIsomorphism.imageGroup(obj));
            else
                c = self.nice.centralizer(self.niceIsomorphism.imageElement(obj));
            end
            c1 = self.niceIsomorphism.preimageGroup(c);
        end

        function res = closure(self, varargin)
            if isempty(varargin)
                res = self;
                return
            end
            rest = cell(1, length(varargin));
            [iso, selfImage, rest{:}] = self.type.niceImages(self, varargin{:});
            nice = selfImage.closure(rest{:});
            res = self.type.groupFromNiceImage(nice, iso);
        end

        function c = complexCharacterTable(self)
            c = self.nice.complexCharacterTable.imap(self.niceIsomorphism.inverse);
        end

        function C = conjugacyClasses(self)
            C = self.cached('conjugacyClasses', @() self.nice.conjugacyClasses.imap(self.niceIsomorphism.inverse));
        end

        function C = conjugacyClass(self, g, varargin)
            args = struct('isCanonical', false, 'centralizer', []);
            args = replab.util.populateStruct(args, varargin);
            g1 = self.niceIsomorphism.imageElement(g);
            if isempty(args.centralizer)
                cc1 = self.nice.conjugacyClass(g1, 'isCanonical', args.isCanonical);
            else
                c1 = args.centralizer.imap(self.niceIsomorphism);
                cc1 = self.nice.conjugacyClass(g1, 'isCanonical', args.isCanonical, 'centralizer', c1);
            end
            C = replab.gen.ConjugacyClass(self.type, cc1, self.niceIsomorphism, self);
        end

        function sub = derivedSubgroup(self)
            sub = self.niceIsomorphism.preimageGroup(self.nice.derivedSubgroup);
        end

        function c = doubleCoset(self, element, rightSubgroup, varargin)
            args = struct('group', [], 'isCanonical', false);
            args = replab.util.populateStruct(args, varargin);
            leftSubgroup = self;
            group = args.group;
            if isempty(group)
                group = leftSubgroup.closure(rightSubgroup, element);
            end
            [iso, group1, leftSubgroup1, element1, rightSubgroup1] = self.type.niceImages(group, leftSubgroup, element, rightSubgroup);
            c1 = leftSubgroup1.doubleCoset(element1, rightSubgroup1, 'group', group1, 'isCanonical', args.isCanonical);
            representative = iso.preimageElement(c1.representative);
            c = replab.gen.DoubleCoset(self.type, c1, iso, representative, leftSubgroup, rightSubgroup, group);
        end

        function c = doubleCosets(self, leftSubgroup, rightSubgroup)
            group = self;
            [iso, group1, leftSubgroup1, rightSubgroup1] = self.type.niceImages(group, leftSubgroup, rightSubgroup);
            c1 = group1.doubleCosets(leftSubgroup1, rightSubgroup1);
            c = replab.gen.DoubleCosets(c1, iso, group, leftSubgroup, rightSubgroup);
        end

        function o = elementOrder(self, g)
            o = self.nice.elementOrder(self.niceIsomorphism.imageElement(g));
        end

        function e = exponent(self)
            e = self.nice.exponent;
        end

        function l = factorizeFlat(self, element)
            if isa(element, 'replab.FiniteSet')
                l = self.nice.factorizeFlat(element.imap(self.niceIsomorphism));
            else
                l = self.nice.factorizeFlat(self.niceIsomorphism.imageElement(element));
            end
        end

        function R = fastRecognize(self)
            R = self.nice.fastRecognize;
            if ~isempty(R)
                R = R.andThen(self.niceIsomorphism.inverse);
            end
        end

        function C = findLeftConjugations(self, s, t, varargin)
            args = struct('sCentralizer', [], 'tCentralizer', []);
            args = replab.util.populateStruct(args, varargin);
            sCentralizer = args.sCentralizer;
            if isempty(sCentralizer)
                sCentralizer = self.centralizer(s);
            end
            if isempty(args.tCentralizer)
                [iso, self1, s1, t1, sCentralizer1] = self.type.niceImages(self, s, t, sCentralizer);
                res1 = self1.findLeftConjugations(s1, t1, 'sCentralizer', sCentralizer1);
            else
                [iso, self1, s1, t1, sCentralizer1, tCentralizer1] = self.type.niceImages(self, s, t, sCentralizer, args.tCentralizer);
                res1 = self1.findLeftConjugations(s1, t1, 'sCentralizer', sCentralizer1, 'tCentralizer', tCentralizer1);
            end
            C = res1.imap(iso.inverse);
        end

        function names = generatorNames(self)
            names = self.nice.generatorNames;
        end

        function sub = intersection(self, rhs)
            [iso, lhs1, rhs1] = self.type.niceImages(self, rhs);
            sub1 = lhs1.intersection(rhs1);
            sub = self.type.groupFromNiceImage(sub1, iso);
        end

        function res = isCyclic(self)
            res = self.nice.isCyclic;
        end

        function res = isSimple(self)
            res = self.nice.isSimple;
        end

        function l = knownComplexCharacterTable(self)
            l = self.nice.knownComplexCharacterTable;
        end

        function res = knownOrder(self)
            res = self.nice.knownOrder;
        end

        function l = knownRealCharacterTable(self)
            l = self.nice.knownRealCharacterTable;
        end

        function res = knownRelators(self)
            res = self.nice.knownRelators;
        end

        function c = leftCoset(self, element, varargin)
            if self.isNormalizedBy(element)
                c = self.normalCoset(element, varargin{:});
                return
            end
            args = struct('group', [], 'isCanonical', false);
            args = replab.util.populateStruct(args, varargin);
            subgroup = self;
            if isempty(args.group)
                group = self.closure(element);
            else
                group = args.group;
            end
            [iso, element1, subgroup1, group1] = self.type.niceImages(element, subgroup, group);
            nice = subgroup1.leftCoset(element1, 'isCanonical', args.isCanonical, 'group', group1);
            representative = iso.preimageElement(nice.representative);
            c = replab.gen.LeftCoset(self.type, nice, iso, representative, subgroup, group);
        end

        function c = leftCosets(self, subgroup)
            group = self;
            [iso, group1, subgroup1] = self.type.niceImages(group, subgroup);
            c1 = group1.leftCosets(subgroup1);
            c = replab.gen.LeftCosets(c1, iso, group, subgroup);
        end

        function res = normalClosure(self, obj)
            [iso, self1, obj1] = self.type.niceImages(self, obj);
            res1 = self1.normalClosure(obj1);
            res = self.type.groupFromNiceImage(res1, iso);
        end

        function c = normalCoset(self, element, varargin)
            assert(self.isNormalizedBy(element));
            args = struct('group', [], 'isCanonical', false);
            args = replab.util.populateStruct(args, varargin);
            subgroup = self;
            if isempty(args.group)
                group = self.closure(element);
            else
                group = args.group;
            end
            [iso, element1, subgroup1, group1] = self.type.niceImages(element, subgroup, group);
            nice = subgroup1.normalCoset(element1, 'isCanonical', args.isCanonical, 'group', group1);
            representative = iso.preimageElement(nice.representative);
            c = replab.gen.NormalCoset(self.type, nice, iso, representative, subgroup, group);
        end

        function c = normalCosets(self, subgroup)
            group = self;
            [iso, group1, subgroup1] = self.type.niceImages(group, subgroup);
            c1 = group1.normalCosets(subgroup1);
            c = replab.gen.NormalCosets(c1, iso, group, subgroup);
        end

        function o = order(self)
            o = self.nice.order;
        end

        function iso = orderPreservingPermutationIsomorphism(self)
            if isa(self.niceIsomorphism.target, 'replab.PermutationGroup')
                iso = self.niceIsomorphism.restrictedSource(self);
            else
                iso = self.niceIsomorphism.restrictedSource(self);
                iso = iso.andThen(self.nice.orderPreservingPermutationIsomorphism);
            end
        end

        function iso = permutationIsomorphism(self)
            if isa(self.niceIsomorphism.target, 'replab.PermutationGroup')
                iso = self.niceIsomorphism.restrictedSource(self);
            else
                iso = self.niceIsomorphism.restrictedSource(self);
                iso = iso.andThen(self.nice.permutationIsomorphism);
            end
        end

        function c = realCharacterTable(self)
            c = self.nice.realCharacterTable.imap(self.niceIsomorphism.inverse);
        end

        function R = relatorsFlat(self)
            R = self.nice.relatorsFlat;
        end

        function c = rightCoset(self, element, varargin)
            if self.isNormalizedBy(element)
                c = self.normalCoset(element, varargin{:});
                return
            end
            args = struct('group', [], 'isCanonical', false);
            args = replab.util.populateStruct(args, varargin);
            subgroup = self;
            if isempty(args.group)
                group = self.closure(element);
            else
                group = args.group;
            end
            [iso, element1, subgroup1, group1] = self.type.niceImages(element, subgroup, group);
            nice = subgroup1.rightCoset(element1, 'isCanonical', args.isCanonical, 'group', group1);
            representative = iso.preimageElement(nice.representative);
            c = replab.gen.RightCoset(self.type, nice, iso, representative, subgroup, group);
        end

        function c = rightCosets(self, subgroup)
            group = self;
            [iso, group1, subgroup1] = self.type.niceImages(group, subgroup);
            c1 = group1.rightCosets(subgroup1);
            c = replab.gen.RightCosets(c1, iso, group, subgroup);
        end

        function iso = regularIsomorphism(self)
            iso = self.niceIsomorphism.restrictedSource(self).andThen(self.nice.regularIsomorphism);
        end

        function rep = regularRep(self)
            rep = self.niceIsomorphism.restrictedSource(self).andThen(self.nice.regularRep);
        end

        function setComplexCharacterTable(self, table)
            self.nice.setComplexCharacterTable(table.imap(self.niceIsomorphism));
        end

        function setConjugacyClasses(self, classes)
            self.nice.setConjugacyClasses(classes.imap(self.niceIsomorphism));
        end

        function setRealCharacterTable(self, table)
            self.nice.setRealCharacterTable(table.imap(self.niceIsomorphism));
        end

    end

    methods % Bugfix for Octave method selection

        function b = contains(self, el)
            b = contains@replab.gen.FiniteSet(self, el);
        end

        function E = elements(self)
            E = elements@replab.gen.FiniteSet(self);
        end

        function E = elementsSequence(self)
            E = elementsSequence@replab.gen.FiniteSet(self);
        end

        function s = nElements(self)
            s = nElements@replab.gen.FiniteSet(self);
        end

        function r = representative(self)
            r = representative@replab.gen.FiniteSet(self);
        end

        function S = setProduct(self)
            S = setProduct@replab.gen.FiniteSet(self);
        end

        function b = isequal(lhs, rhs)
            b = isequal@replab.FiniteGroup(lhs, rhs);
        end

    end

    methods

        % Workaround for Octave bug

        function res = eq(self, rhs)
            res = replab.finite.equality(self, rhs);
        end

        function res = ne(self, rhs)
            res = ~replab.finite.equality(self, rhs);
        end

    end

end
