classdef TypeToPermIsomorphism < replab.FiniteIsomorphism
% A morphism from a subgroup of a finite group type to a permutation group

    methods

        function t = sourceType(self)
        % Returns the type of the finite group which is the source of this isomorphism
            t = self.source.type;
        end
        function l = sourceContains(self, s)
        % Returns whether the source of this isomorphism contains the given type element
        %
        % Args:
        %   s (element of `.sourceType`): Element to check
        %
        % Returns:
        %   logical: True if source contains the given element
            error('Abstract');
        end

    end

end
