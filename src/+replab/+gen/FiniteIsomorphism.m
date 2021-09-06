classdef FiniteIsomorphism < replab.FiniteIsomorphism
% An isomorphism from a subgroup of a finite group type to a finite group (Abstract)

    methods

        function self = FiniteIsomorphism(sourceFun, targetType)
        % Constructs a finite isomorphism
        %
        % Args:
        %   sourceFun (function_handle): Function handle that takes this partially constructed object as argument and returns the source group
        %   targetType (`+replab.FiniteGroupType`): Type of target finite group
            source = sourceFun(self);
            self.source = source;
            targetGenerators = cellfun(@(g) self.imageElement(g), source.generators, 'uniform', 0);
            self.target = targetType.groupWithGenerators(targetGenerators);

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
