classdef FiniteGroup < replab.FiniteGroup & replab.gen.FiniteSet
% A generic finite group

    methods

        function self = FiniteGroup(type, nice, niceIsomorphism, generators)
        % Constructs a nice finite group
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
            self.identity = type.identity;
            self.representative_ = type.identity;
            self.generators = generators;
        end

    end

% $$$     methods (Access = protected) % Implementations
% $$$
% $$$         function C = computeConjugacyClasses(self)
% $$$             nc = self.genericGroup.conjugacyClasses;
% $$$             C = nc.imap(self.genericIsomorphism.inverse);
% $$$         end
% $$$
% $$$         function sub = computeDerivedSubgroup(self)
% $$$             sub = self.genericIsomorphism.preimageGroup(self.genericGroup.derivedSubgroup);
% $$$         end
% $$$
% $$$         function E = computeElements(self)
% $$$             atFun = @(ind) self.genericIsomorphism.preimageElement(self.genericGroup.elements.at(ind));
% $$$             findFun = @(el) self.genericGroup.elements.find(self.genericIsomorphism.imageElement(el));
% $$$             E = replab.IndexedFamily.lambda(self.order, atFun, findFun);
% $$$         end
% $$$
% $$$         function res = computeIsCyclic(self)
% $$$             res = self.genericGroup.isCyclic;
% $$$         end
% $$$
% $$$         function res = computeIsSimple(self)
% $$$             res = self.genericGroup.isSimple;
% $$$         end
% $$$
% $$$         function order = computeOrder(self)
% $$$             order = self.genericGroup.order;
% $$$         end
% $$$
% $$$         function dec = computeSetProduct(self)
% $$$             prmD = self.genericGroup.setProduct;
% $$$             sets1 = cellfun(@(T) cellfun(@(t) self.genericIsomorphism.preimageElement(t), T, 'uniform', 0), prmD.sets, 'uniform', 0);
% $$$             dec = replab.SetProduct(self, sets1, true);
% $$$         end
% $$$
% $$$     end

    methods % Implementation

        % Domain

        function b = eqv(self, x, y)
            b = self.type.eqv(x, y);
        end

        function g = sample(self)
            g = self.niceIsomorphism.preimageElement(self.nice.sample);
        end

        % FiniteSet

        % Monoid

        function z = compose(self, x, y)
            z = self.type.compose(x, y);
        end

        % Group

        function y = inverse(self, x)
            y = self.type.inverse(x);
        end

        % FiniteGroup

% $$$         function A = abstractGroup(self, generatorNames)
% $$$             if nargin < 2 || isempty(generatorNames)
% $$$                 generatorNames = self.generatorNames;
% $$$             end
% $$$             A = self.genericGroup.abstractGroup(generatorNames);
% $$$         end
% $$$
% $$$         function m = abstractMorphism(self, generatorNames)
% $$$             if nargin < 2 || isempty(generatorNames)
% $$$                 generatorNames = self.generatorNames;
% $$$             end
% $$$             m = self.genericIsomorphism.andThen(self.genericGroup.abstractMorphism(generatorNames));
% $$$         end

        % Group elements

        function b = contains(self, g)
            if self.niceIsomorphism.sourceContains(g)
                b = self.nice.contains(self.niceIsomorphism.imageElement(g));
            else
                b = false;
            end
        end

        function o = elementOrder(self, g)
            o = self.nice.elementOrder(self.niceIsomorphism.imageElement(g));
        end

        function l = factorizeLetters(self, element)
            l = self.nice.factorizeLetters(self.niceIsomorphism.imageElement(element));
        end

        % Group properties

        function a = abelianInvariants(self)
            a = self.nice.abelianInvariants;
        end

        function e = exponent(self)
            e = self.nice.exponent;
        end

        function o = order(self)
            o = self.nice.order;
        end

        % Construction of groups

% $$$         function res1 = closure(self, obj)
% $$$             if isa(obj, 'replab.FiniteGroup')
% $$$                 % if one group contains the other
% $$$                 if self.isSubgroupOf(obj)
% $$$                     res1 = obj;
% $$$                     return
% $$$                 end
% $$$                 if obj.isSubgroup(self)
% $$$                     res1 = self;
% $$$                     return
% $$$                 end
% $$$                 % otherwise do the computation
% $$$                 res = self.genericGroup.closure(self.type.niceMorphism.imageGroup(obj));
% $$$             else
% $$$                 % if the group already contains the element
% $$$                 if self.contains(obj)
% $$$                     res1 = self;
% $$$                     return
% $$$                 end
% $$$                 % otherwise do the computation
% $$$                 res = self.niceGroup.closure(self.type.niceImage(obj));
% $$$             end
% $$$             res1 = self.type.niceMorphism.preimageGroup(res);
% $$$         end
% $$$
% $$$         function res1 = normalClosure(self, obj)
% $$$             if isa(obj, 'replab.FiniteGroup')
% $$$                 res = self.niceGroup.normalClosure(self.type.niceMorphism.imageGroup(obj));
% $$$             else
% $$$                 res = self.niceGroup.normalClosure(self.type.niceImage(obj));
% $$$             end
% $$$             res1 = self.type.niceMorphism.preimageGroup(res);
% $$$         end

        % Subgroups

        function sub = subgroupWithGenerators(self, generators, order)
            if nargin < 3
                order = [];
            end
            if length(generators) == self.nGenerators && all(arrayfun(@(i) self.eqv(self.generator(i), generators{i}), 1:length(generators)))
                sub = self;
            else
                niceGroup = self.niceGroup.subgroupWithGenerators(cellfun(@(g) self.niceImage(g), generators, 'uniform', 0), order);
                sub = self.niceSubgroup(generators, order, niceGroup);
            end
        end

        function sub1 = centralizer(self, obj)
            if isa(obj, 'replab.FiniteGroup')
                sub = self.niceGroup.centralizer(self.niceMorphism.imageGroup(obj));
            else
                sub = self.niceGroup.centralizer(self.niceImage(obj));
            end
            sub1 = self.niceMorphism.preimageGroup(sub);
        end

        function sub1 = intersection(self, rhs)
            sub1 = self.niceMorphism.preimageGroup(self.niceGroup.intersection(self.niceMorphism.imageGroup(rhs)));
        end

        % Cosets

        function B = findLeftConjugations(self, s, t, sCentralizer, tCentralizer)
            s = self.niceImage(s);
            t = self.niceImage(t);
            if nargin < 4 || isempty(sCentralizer)
                sCentralizer = [];
            else
                sCentralizer = self.niceMorphism.imageGroup(sCentralizer);
            end
            if nargin < 5 || isempty(tCentralizer)
                tCentralizer = [];
            else
                tCentralizer = self.niceMorphism.imageGroup(tCentralizer);
            end
            B = self.niceGroup.findLeftConjugations(s, t, sCentralizer, tCentralizer);
            group = self.niceMorphism.preimageGroup(B.group);
            canRep = self.niceMorphism.preimageElement(B.representative);
            B = replab.LeftCoset(group, canRep, self);
        end

        % Relation to other groups

        % Representations

        function rep = regularRep(self)
            rep = self.niceMorphism.andThen(self.niceGroup.regularRep);
        end

    end

    methods (Access = protected)

        % Morphisms

        function m = morphismByImages_(self, target, preimages, images, imageElementFun)
            first = self.niceMorphism; % maps this to the perm group
            preimagesNG = cellfun(@(g) self.niceMorphism.imageElement(g), preimages, 'uniform', 0);
            second = self.niceGroup.morphismByImages(target, 'preimages', preimagesNG, 'images', images, 'nChecks', 0); % from the perm group to the images
            if isa(second, 'replab.FiniteMorphism')
                m = replab.mrp.FiniteComposition(second, first, imageElementFun);
            else
                m = replab.mrp.Composition(second, first, imageElementFun);
            end
        end

    end

end
