classdef PlainAutomorphism < replab.FiniteIsomorphism
% Describes an automorphism in a finite group through an automorphism of its permutation realization

    properties (SetAccess = protected)
        object % (`+replab.FiniteGroup`): Group
        generatorImages % (cell(1,\*) of elements of `.object`): Images of generators defining an automorphism
    end

    methods

        function self = PlainAutomorphism(object, generatorImages)
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
            m = replab.nfg.PlainAutomorphism(self.object, generatorImages1);
        end

    end

    methods % Implementations

        % Str

        function s = shortStr(self, maxColumns)
            images = arrayfun(@(i) [replab.shortStr(self.object.generator(i)) ' -> ' replab.shortStr(self.generatorImages{i})], 1:self.object.nGenerators, 'uniform', 0);
            s = ['PlainAutomorphism(' strjoin(images, ', ') ')'];
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
        % Creates an automorphism of a group by conjugation
        %
        % Args:
        %   h (Group element): Element of the same type as ``object.type``
        %   object (`+replab.FiniteGroup`): Finite group being conjugated
        %
        % Returns:
        %   `.PlainAutomorphism`: A conjugation automorphism
            generatorImages = cellfun(@(g) object.type.leftConjugate(h, g), object.generators, 'uniform', 0);
            assert(all(cellfun(@(g) object.contains(g), generatorImages)), 'The conjugation does not define an automorphism');
            a = replab.nfg.PlainAutomorphism(object, generatorImages);
        end

        function a = identity(object)
        % Returns the identity automorphism
        %
        % Args:
        %   object (`+replab.FiniteGroup`): Finite group
        %
        % Returns:
        %   `.PlainAutomorphism`: Identity automorphism
            generatorImages = object.generators;
            a = replab.nfg.PlainAutomorphism(object, generatorImages);
        end

    end

end
