classdef FiniteMorphism < replab.Morphism
% Describes a morphism between finite groups

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.FiniteMorphismLaws(self);
        end

    end

    methods (Access = protected)

        function K = computeKernel(self)
            K = self.preimageGroup(self.target.trivialSubgroup);
        end

        function I = computeImage(self)
            I = self.imageGroup(self.source);
        end

        function I = computeImageSourceGenerators(self)
            I = cellfun(@(s) self.imageElement(s), self.source.generators, 'uniform', 0);
        end

    end

    methods % Conversions

        function m = restrictedSource(self, s1)
        % Returns a finite morphism with its source restricted
        %
        % Note that the return type may be more precise.
        %
        % Args:
        %   newSource (`.FiniteGroup`): Subgroup of `.source`
        %
        % Returns:
        %   `.FiniteMorphism`: The restricted morphism
            assert(newSource.isSubgroupOf(self.source), 'Must restrict to a subgroup of source');
            if self.source.isSubgroupOf(newSource)
                m = self; % we have newSource == self.source
            else
                images = cellfun(@(g) self.imageElement(g), s1.generators, 'uniform', 0);
                m = s1.morphismByImages(self.target, images);
            end
        end

        function m = toIsomorphism(self)
        % If this morphism is injective, returns an isomorphism such that the target is restricted to the image
        %
        % Returns:
        %   `.FiniteIsomorphism`: The isomorphism
            assert(self.kernel.isTrivial);
            m = replab.fm.FiniteIsomorphismWrapper(self);
        end

    end

    methods % Preimages

        function K = kernel(self)
        % Returns the kernel of this morphism
        %
        % Returns:
        %   `+replab.FiniteGroup`: Maximal subgroup of `.source` with trivial image
            K = self.cached('kernel', @() self.computeKernel);
        end

        function s = preimageRepresentative(self, t)
        % Returns an arbitrary preimage of the given element
        %
        % Returns an ``s`` such that ``self.imageElement(s) == t`` .

        % Args:
        %   t (element of `.target`): Element to compute the preimage of
        %
        % Returns:
        %   element of `.source`: Preimage representative
            error('Abstract');
        end

        function S = preimagesElement(self, t)
        % Returns the set of all source elements that map to a given element
        %
        % Args:
        %   t (element of `.target`): Element to compute the preimages of
        %
        % Returns:
        %   `.FiniteSet`: Set of source elements
            S = self.source.normalCoset(self.kernel, self.preimageRepresentative(t));
        end

        function S = preimageGroup(self, T)
        % Returns the group of source elements that map to a given group of target elements
        %
        % Args:
        %   T (`.FiniteGroup`): Subgroup of `.target`
        %
        % Returns:
        %   `.FiniteGroup`: Subgroup of `.source`
            error('Abstract');
        end

    end

    methods % Images

        function I = imageSourceGenerators(self)
        % Returns the images of the source generators
        %
        % Returns:
        %   cell(1,\*) of elements of target: Generator images
            I = self.cached('imageSourceGenerators', @() self.computeImageSourceGenerators);
        end

        function I = image(self)
        % Returns the image of this morphism
        %
        % Returns:
        %   `.FiniteGroup`: Subgroup of `.target` generated by the image of `.source`
            I = self.cached('image', @() self.computeImage);
        end

        function T = imageGroup(self, S)
        % Computes the image of a group
        %
        % Args:
        %   S (`.FiniteGroup`): Group to compute the image of, subgroup of `.source`
        %
        % Returns:
        %   `.FiniteGroup`: Subgroup of `.target`
            images = cellfun(@(g) self.imageElement(g), S.generators, 'uniform', 0);
            T = self.target.subgroup(images);
        end

    end

end
