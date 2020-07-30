classdef FiniteIsomorphism < replab.Isomorphism & replab.FiniteMorphism
% Describes an isomorphism between finite groups

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.FiniteIsomorphismLaws(self);
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

    methods % Implementations

        % FiniteMorphism

        function s = preimageRepresentative(self, t)
            s = t;
        end

        function S = preimageGroup(T)
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
