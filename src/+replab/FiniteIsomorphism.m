classdef FiniteIsomorphism < replab.Isomorphism & replab.FiniteMorphism
% Describes an isomorphism between finite groups

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.FiniteIsomorphismLaws(self);
        end

    end

    methods (Access = protected)

        function I = computeInverse(self)
            I = replab.fm.FiniteInverse(self);
        end

        function K = computeKernel(self)
            K = self.source.trivialSubgroup;
        end

    end

    methods

        function S = preimagesElement(self, t)
            S = self.normalCoset(self.source, self.kernel, self.preimageElement(t));
        end

        function T = imageGroup(self, S)
            images = cellfun(@(g) self.imageElement(g), S.generators, 'uniform', 0);
            T = self.target.subgroupWithGenerators(images); % upgraded
        end

        function S = preimageGroup(T)
            S = self.inverse.imageGroup(T);
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
            m = replab.fm.Identity(group);
        end

    end

end
