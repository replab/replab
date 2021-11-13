classdef StaticFiniteGroupType < replab.gen.FiniteGroupType
% A static finite group type uses the same isomorphism for every finite group it creates
%
% The adjective "static", in this context, does not have anything to do with the static/dynamic typing
% terminology used in computer programming.
%
% To ease the implementation such static finite group types, we put the isomorphim abstract methods
% to compute preimages and images in this class, such that the user only needs to subclass
% `+replab.+gen.StaticFiniteGroupType` to define a new group type in RepLAB.
%
% There is a single challenge

    properties (SetAccess = protected)
        niceType % (`+replab.FiniteGroupType`): Nice type
        isomorphism % (`.NiceIsomorphism`): Isomorphism whose source contains all elements of this type
        parentGroup % (`+replab.FiniteGroup`): Group that contains all the type elements, equal to ``self.isomorphism.source``
    end

    methods (Access = protected)

        function finishConstruction(self, sourceGenerators, sourceArgs, niceType)
        % Finish the construction of this static group type
        %
        % The code below should be called at the end of the subclass constructor.
        %
        % It constructs the nice isomorphism used for all groups of this type, and does
        % the following:
        %
        % * It translates the given source generators to construct the
        %   nice group / isomorphism target.
        %
        % * It constructs the default nice isomorphism implementation, which will itself
        %
        % At the call site, the finite group type object will be partially constructed.
        % This method needs that `.imageElement` and `.makeIsomorphismSourceGroup` work
        % at that point.
        %
        % This code path is a bit convoluted, as required by the construction of a graph of
        % objects with circular references:
        %
        % * This group type references both `.parentGroup` and `.isomorphism`.
        %
        % * The group `.parentGroup` has a `+replab.+gen.FiniteGroup.niceIsomorphism` property
        %   equal to `.isomorphism`, and a `+replab.+gen.FiniteGroup.nice` property equal to
        %   `+replab.+gen.NiceIsomorphism.target`.
        %
        % * This `.isomorphism` property `+replab.+gen.NiceIsomorphism.source` is equal to
        %   this type `.parentGroup` property.
        %
        % Args:
        %   sourceGenerators (cell(1,\*) of elements of this type): Generators of the parent group
        %   sourceArgs (cell(1,\*) of key/value pairs): Contains generator names, group order, relators if known
        %   niceType (`+replab.FiniteGroupType`): Type of the group where computations are delegated
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

        function G = makeParentGroup(self, generators, nice, niceIsomorphism)
        % Creates the source group of the nice isomorphism (internal)
        %
        % Has the same calling convention than `.makeGenericGroup` but can be overriden
            G = self.makeGenericGroup(generators, nice, niceIsomorphism);
        end

        function setParentGroup(self, parentGroup)
        % Sets the parent group (internal)
        %
        % Hack used during construction of `+replab.AbstractGroup`, to sovle a problem with
        % circular references.
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

        function iso = constructNiceIsomorphism(self, elements)
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
