classdef FiniteAutomorphismGroup < replab.NiceFiniteGroup
% Describes the automorphism group of a group
%
% We store the semidirect product of the outer automorphism group with the inner automorphism group.
%
% The outer automorphism group is described by an isomorphic permutation group, with a map (`outerAutomorphisms`) that provides normal
% coset representative images.
%
% The inner automorphism group is described by coset representatives of the group we study by its center.
%
% This code has bugs when the outer automorphism group is not trivial!
%
% Example:
%   >>> S = replab.S(6);
%   >>> outerAutomorphisms = replab.S(2).morphismByImages(replab.nfg.PlainAutomorphismGroup(S), {S.morphismByImages(S, {[6 1 5 4 3 2], [2 1 4 3 6 5]})});
%   >>> A = replab.FiniteAutomorphismGroup(S, outerAutomorphisms);
%
    properties (Access = protected)
        innerAutomorphisms % (`.NormalCosets`): Cosets of `.object` over ``object.center``
        outerAutomorphisms % (`.FiniteMorphism`): Homomorphism of a permutation group whose image is the outer automorphism group
    end

    properties (SetAccess = protected)
        object % (`.FiniteGroup`): Group
    end

    methods

        function self = FiniteAutomorphismGroup(object, outerAutomorphisms)
            self.object = object;
            self.innerAutomorphisms = object.normalCosetsOf(object.center);
            self.outerAutomorphisms = outerAutomorphisms;
            self.identity = replab.FiniteAutomorphism(self, outerAutomorphisms.source.identity, object.identity);
            mask = cellfun(@(g) object.center.contains(g), object.generators);
            innerGenerators = object.generators(mask);
            outerGenerators = outerAutomorphisms.source.generators;
            innerGenerators = cellfun(@(ig) replab.FiniteAutomorphism(self, outerAutomorphisms.source.identity, ig), innerGenerators, 'uniform', 0);
            outerGenerators = cellfun(@(og) replab.FiniteAutomorphism(self, og, self.object.identity), outerGenerators, 'uniform', 0);
            self.generators = horzcat(outerGenerators, innerGenerators);
            self.type = self;
        end

        function s = objectSet(self)
            s = self.cached('objectSet', @() self.computeObjectSet);
        end

        function o = order(self)
            o = self.innerAutomorphisms.nElements * self.outerAutomorphisms.source.order;
        end

        function e = elements(self)
            e = self.cache('elements', @() self.computeElements);
        end

        function g1 = outerAction(self, outerPreimage, g, applyInverse)
            if nargin < 4
                applyInverse = false;
            end
            m = self.outerAutomorphisms.imageElement(outerPreimage);
            if applyInverse
                g2 = m.preimageElement(g);
            else
                g2 = m.imageElement(g);
            end
            ab = self.outerAutomorphisms.source.abstractGroupIsomorphism;
            w = ab.target.toLetters(ab.imageElement(outerPreimage));
            if applyInverse
                w = -fliplr(w);
            end
            g1 = g;
            imgs = self.outerAutomorphisms.imageSourceGenerators;
            for i = length(w):-1:1 % left action convention
                wi = w(i);
                img = imgs{abs(wi)};
                if wi > 0
                    g1 = img.imageElement(g1);
                else
                    g1 = img.preimageRepresentative(g1);
                end
            end
            g1 = self.innerAutomorphisms.cosetRepresentative(g1);
            g2 = self.innerAutomorphisms.cosetRepresentative(g2);
            assert(self.object.eqv(g1, g2));
        end

    end

    methods (Access = protected)

        function s = computeObjectSet(self)
            permGrp = self.object.niceMorphism.image;
            s = replab.perm.Set.fromPermutationGroup(permGrp);
            s.sort;
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = self.outerAutomorphisms.source.eqv(x.outerPreimage, y.outerPreimage) && self.object.eqv(x.innerRepresentative, y.innerRepresentative);
        end

        function s = sample(self, x, y)
            s = replab.FiniteAutomorphism(self, self.outerAutomorphisms.source.sample, self.innerAutomorphisms.cosetRepresentative(self.object.sample));
        end

        % Monoid

        function x3 = compose(self, x1, x2)
            o1 = x1.outerPreimage;
            i1 = x1.innerRepresentative;
            o2 = x2.outerPreimage;
            i2 = x2.innerRepresentative;
            % g -> o3(i3 * g * i3^-1)
            % g -> o1(o2(i3 * g * i3^-1))
            % g -> o1(i1 * o2(i2 * g * i2^-1) * i1^-1)
            % g -> o1(o2(o2^-1(i1) * i2 * g * i2^-1 * o2^-1(i1)^-1))
            % o3 = g -> o1(o2(g))
            % i3 = o2^-1(i1) * i2
            o3 = self.outerAutomorphisms.source.compose(o1, o2);
            i3 = self.innerAutomorphisms.cosetRepresentative(self.object.compose(self.outerAction(o2, i1, true), i2));
            x3 = replab.FiniteAutomorphism(self, o3, i3);
            for i = 1:self.object.nGenerators
                g = self.object.generator(i);
                img1 = x1.imageElement(x2.imageElement(g));
                img2 = x3.imageElement(g);
                assert(self.object.eqv(img1, img2));
            end
        end

        function x2 = inverse(self, x1)
            o1 = x1.outerPreimage;
            i1 = x1.innerRepresentative;
            o2 = self.outerAutomorphisms.source.inverse(o1);
            i2 = self.innerAutomorphisms.cosetRepresentative(self.object.inverse(self.outerAction(o2, i1, true)));
            x2 = replab.FiniteAutomorphism(self, o2, i2);
        end

        % NiceFiniteGroup

        function img = niceImage(self, g)
            n = double(self.object.order);
            img = zeros(1, n);
            for i = 1:n
                img(i) = self.object.elements.find(g.imageElement(self.object.elements.at(i)));
            end
        end

    end

end
