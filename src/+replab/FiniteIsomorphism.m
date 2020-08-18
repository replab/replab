classdef FiniteIsomorphism < replab.Isomorphism & replab.FiniteMorphism
% Describes an isomorphism between finite groups
%
% Adds the guarantee that `.target`/`.image` has for generators the images of the generators of `.source`

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.FiniteIsomorphismLaws(self);
        end

    end

    methods (Access = protected)

        function I = computeImage(self)
            I = self.target; % due to self being an isomorphism
        end

        function I = computeInverse(self)
            I = replab.fm.FiniteInverse(self);
        end

        function K = computeKernel(self)
            K = self.source.trivialSubgroup;
        end

    end

    methods % Implementations

        % FiniteMorphism

        function s = preimageRepresentative(self, t)
            s = self.preimageElement(t);
        end

        function S = preimageGroup(self, T)
            S = self.inverse.imageGroup(T);
        end

        function T = imageGroup(self, S)
            images = cellfun(@(g) self.imageElement(g), S.generators, 'uniform', 0);
            T = self.target.subgroupWithGenerators(images); % do not need to check for non-generators
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
            m = replab.fm.FiniteIdentity(group);
        end

    end

end
