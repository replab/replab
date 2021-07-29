classdef GenericIsomorphism < replab.FiniteIsomorphism
% An isomorphism from a subgroup of a finite group type to a finite group

    methods

        function self = GenericIsomorphism(source, target)
            self.source = source;
            self.target = target;
            self.torusMap = [];
        end

        function t = sourceType(self)
        % Returns the type of the finite group which is the source of this morphism
            t = self.source.type;
        end

        function l = sourceContains(self, s)
        % Returns whether the source of this morphism contains the given type element
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
