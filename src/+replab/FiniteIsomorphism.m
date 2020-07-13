classdef FiniteIsomorphism < replab.FiniteMorphism
% Describes an isomorphism between finite groups

    methods

        function K = computeKernel(self)
            K = self.source.trivialSubgroup;
        end

        %function S = preimagesElement(self, t)
        %    S = {self.preimageElement(t)};
        %end

        function T = imageGroup(self, S)
            images = cellfun(@(g) self.imageElement(g), S.generators, 'uniform', 0);
            T = self.target.subgroupWithGenerators(images); % upgraded
        end

        function s = preimageElement(t)
            error('Abstract');
        end

        function S = preimageGroup(T)
            S = self.inverse.imageGroup(T);
        end

        function I = inverse(self)
            I = replab.fm.Inverse(self);
        end

        function res = compose(self, applyFirst)
            if isa(applyFirst, 'replab.FiniteIsomorphism')
                m = replab.fm.IsoComposition(self, applyFirst);
            else
                m = replab.fm.Composition(self, applyFirst);
            end
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
