classdef FiniteIsomorphism < replab.Isomorphism & replab.FiniteMorphism
% Describes an isomorphism between finite groups

    methods

        function l = preservesTypeOrder(self)
        % Returns whether, for sure, this isomorphism preserves element ordering
        %
        % * When this method returns true, we must have, for ``x`` and ``y`` in `.source`,
        %   ``self.source.type.compare(x, y) == self.target.type.compare(self.imageElement(x), self.imageElement(y))``
        % * When this method returns false, it does not necessarily means that there exists a pair ``(x, y)``
        %   that violates the condition above. It simply means the order preservation is not known.
        %
        % Returns:
        %   logical: Whether, for sure, this isomorphism preserves element ordering
            l = false; % by default, if unknown
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.FiniteIsomorphismLaws(self);
        end

        % FiniteMorphism

        function T = imageGroup(self, S)
            images = cellfun(@(g) self.imageElement(g), S.generators, 'uniform', 0);
            T = self.target.subgroupWithGenerators(images); % do not need to check for non-generators
        end

        function S = preimageGroup(self, T)
            preimages = cellfun(@(g) self.preimageElement(g), T.generators, 'uniform', 0);
            S = self.source.subgroupWithGenerators(preimages); % do not need to check for non-generators
        end

        function s = preimageRepresentative(self, t)
            s = self.preimageElement(t);
        end

        function m = restrictedSource(self, newSource)
            m = replab.mrp.SourceRestrictedFiniteIsomorphism(self, newSource);
        end

        % Str

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteMorphism(self);
        end

    end

    methods (Access = protected)

        function I = computeImage(self)
            I = self.target; % due to self being an isomorphism
        end

        function I = computeInverse(self)
            I = replab.mrp.FiniteInverse(self);
        end

        function K = computeKernel(self)
            K = self.source.trivialSubgroup;
        end

    end

    methods (Static)

        function m = identity(group)
        % Returns the identity morphism from a finite group to itself
        %
        % Args:
        %   group (`.FiniteGroup`): Group
        %
        % Returns:
        %   `.FiniteIsomorphism`: The identity automorphism on the given group
            m = replab.mrp.FiniteIdentity(group);
        end
        
    end

end
