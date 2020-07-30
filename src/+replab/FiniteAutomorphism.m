classdef FiniteAutomorphism < replab.FiniteIsomorphism

    properties (SetAccess = protected)
        parent % (`.FiniteAutomorphismGroup`): Automorphism group this is part of
        outerPreimage % (permutation): Preimage of the outer automorphism part in ``group.outerAutomorphisms``
        innerRepresentative % (element): Representative of a normal coset in ``group.innerAutomorphisms``
    end

    methods

        function self = FiniteAutomorphism(parent, outerPreimage, innerRepresentative)
            self.parent = parent;
            self.source = parent.object;
            self.target = parent.object;
            self.outerPreimage = outerPreimage;
            self.innerRepresentative = innerRepresentative;
        end

    end

    methods % Implementations

        % Str

        function s = shortStr(self, maxColumns)
            images = arrayfun(@(i) [replab.shortStr(self.parent.object.generator(i)) ' -> ' replab.shortStr(self.imageElement(self.parent.object.generator(i)))], 1:self.parent.object.nGenerators, 'uniform', 0);
            s = ['FiniteAutomorphism(' strjoin(images, ', ') ')'];
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.FiniteIsomorphism(self);
            for i = 1:self.parent.object.nGenerators
                g = self.parent.object.generator(i);
                names{1,end+1} = sprintf('imageElement(%s)', replab.shortStr(g));
                values{1, end+1} = self.imageElement(g);
            end
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.FiniteIsomorphism(self);
            names{1, end+1} = 'source';
            names{1, end+1} = 'target';
        end

        % Morphism

        function g2 = imageElement(self, g)
            g1 = self.parent.object.leftConjugate(self.innerRepresentative, g);
            g2 = self.parent.outerAction(self.outerPreimage, g1, false);
        end

        % Isomorphism

        function g2 = preimageElement(self, g)
            g1 = self.parent.outerAction(self.outerPreimage, g, true);
            g2 = self.parent.object.leftConjugate(self.parent.object.inverse(self.innerRepresentative), g1);
        end

    end

end
