classdef StaticFiniteGroupType < replab.gen.FiniteGroupType
% A static finite group type is a type whose number of elements is reasonably bounded
%
% Using the bounded cardinality, we can construct a single isomorphism that is used to create all generic groups.

    properties (SetAccess = protected)
        niceType % (`+replab.FiniteGroupType`): Nice type
        isomorphism % (`.NiceIsomorphism`): Isomorphism whose source contains all elements of this type
        parentGroup % (`+replab.FiniteGroup`): Group that contains all the type elements, equal to ``self.isomorphism.source``
    end

    methods (Access = protected)

        function finishConstruction(self, sourceGenerators, sourceArgs, niceType)
            targetGenerators = cellfun(@(g) self.imageElement(g), sourceGenerators, 'uniform', 0);
            % Verify that we can bump the source properties to the target group
            assert(all(ismember(sourceArgs(1:2:end), {'generatorNames', 'order', 'relators'})));
            targetArgs = sourceArgs;
            % Construct the nice group
            target = niceType.groupWithGenerators(targetGenerators, targetArgs{:});
            self.niceType = niceType;
            self.isomorphism = replab.gen.StaticNiceIsomorphism(self, sourceGenerators, target);
            self.parentGroup = self.isomorphism.source;
        end

    end

    methods

        function G = makeIsomorphismSourceGroup(self, generators, nice, niceIsomorphism)
        % Creates the source group of the nice isomorphism (internal)
        %
        % Has the same calling convention than `.makeGenericGroup` but can be overriden
            G = self.makeGenericGroup(generators, nice, niceIsomorphism);
        end

        function setParentGroup(self, parentGroup)
        % Sets the parent group (internal)
        %
        % Hack used during construction of `+replab.AbstractGroup`, as the group is constructed
        % before the group type.
        %
        % Args:
        %   parentGroup (`+replab.FiniteGroup`): Parent group to set
            self.parentGroup = parentGroup;
        end

    end

    methods % Implementations

        % TotalOrder

        function c = compare(self, x, y)
            c = self.niceType.compare(self.imageElement(x), self.imageElement(y));
        end

        % FiniteGroupType

        function iso = constructIsomorphism(self, elements)
            iso = self.isomorphism;
        end

    end

    methods % Internal methods

        % To be overriden by finite group type implementations

        function t = imageElement(self, s)
        % Computes the (nice) image of an element of this type (internal)
        %
        % Args:
        %   s (element of ``self``): Element of this type
        %
        % Returns:
        %   element of `.niceType`: Image of the element
            error('Abstract');
        end

        function s = preimageElement(self, t)
        % Computes the preimage of an element (internal)
        %
        % Args:
        %   t (element of `.isomorphism` target): Nice image of an element
        %
        % Returns:
        %   element of ``self``: Preimage of the element
            error('Abstract');
        end

    end

end
