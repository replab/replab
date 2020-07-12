classdef ConjugacyClass < replab.Str

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group containing this conjugacy class
        representative % (group element): Representative element of the conjugacy class
        size % (integer): Size of the conjugacy class
    end

    properties (Access = protected)
        morphism % (`.FiniteMorphism`): Morphism that maps permutations back to elements of the original group
        elementImages % (double(\*,\*)): Matrix containing images of the conjugacy class elements
                      %
                      %                  Each column is a permutation representing the image
                      %                  of an element of the conjugacy class through
                      %                  `+replab.NiceFiniteGroup.niceMonomorphismImage`
    end

    methods (Access = protected)

        function self = ConjugacyClass(group, morphism, elementImages)
            self.group = group;
            self.representative = morphism.image(elementImages(:,1)');
            self.morphism = morphism;
            self.elementImages = elementImages;
            self.size = size(self.elementImages, 2);
        end

    end

    methods

        function e = elements(self)
            e = arrayfun(@(i) self.morphism.image(self.elementImages(:,i)'), 1:self.size, 'uniform', 0);
        end

    end

    methods (Static)

        function c = computeAll(group)
            classes = replab.nfg.conjugacyClassesByOrbits(group);
            c = cellfun(@(cl) replab.ConjugacyClass(group, cl), classes, 'uniform', 0);
        end

    end

end
