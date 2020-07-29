classdef Automorphism < replab.FiniteIsomorphism

    properties (SetAccess = protected)
        object % (`.FiniteGroup`): Group
        generatorImages % (cell(1,\*) of elements of `.object`): Images of generators defining an automorphism
    end

    methods

        function self = Automorphism(object, generatorImages)
            self.source = object;
            self.target = object;
            self.object = object;
            self.generatorImages = generatorImages;
        end

        function m = niceAutomorphism(self)
            m = self.cached('niceAutomorphism', @() self.computeNiceAutomorphism);
        end

    end

    methods (Access = protected)

        function m = computeNiceAutomorphism(self)
            permGrp = self.object.niceMorphism.image;
            permImages = cellfun(@(g) self.object.niceMorphism.imageElement(g), self.generatorImages, 'uniform', 0);
            m = replab.fm.PermToPerm(permGrp, permGrp, permImages);
        end

        function m = computeInverse(self)
            generatorImages1 = cellfun(@(g) self.preimageElement(g), self.object.generators, 'uniform', 0);
            m = replab.Automorphism(self.object, generatorImages1);
        end

    end

    methods % Implementations

        % Str

        function s = shortStr(self)
            images = arrayfun(@(i) [replab.shortStr(self.object.generator(i)) ' -> ' replab.shortStr(self.generatorImages{i})], 1:self.object.nGenerators, 'uniform', 0);
            s = ['Automorphism(' strjoin(images, ', ') ')'];
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.FiniteIsomorphism(self);
            for i = 1:self.object.nGenerators
                names{1,end+1} = sprintf('imageElement(%s)', replab.shortStr(self.object.generator(i)));
                values{1, end+1} = self.generatorImages{i};
            end
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.FiniteIsomorphism(self);
            names{1, end+1} = 'source';
            names{1, end+1} = 'target';
            names{1, end+1} = 'generatorImages';
        end

        % Morphism

        function g1 = imageElement(self, g)
            g1 = self.object.niceMorphism.preimageElement(self.niceAutomorphism.imageElement(self.object.niceMorphism.imageElement(g)));
        end

        % Isomorphism

        function g1 = preimageElement(self, g)
            g1 = self.object.niceMorphism.preimageElement(self.niceAutomorphism.preimageRepresentative(self.object.niceMorphism.imageElement(g)));
        end

    end

    methods (Static)

        function a = byConjugation(h, object)
            generatorImages = cellfun(@(g) object.leftConjugate(h, g), object.generators, 'uniform', 0);
            a = replab.Automorphism(object, generatorImages);
        end

        function a = identity(object)
            generatorImages = object.generators;
            a = replab.Automorphism(object, generatorImages);
        end

    end

end
