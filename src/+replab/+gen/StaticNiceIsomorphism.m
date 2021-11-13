classdef StaticNiceIsomorphism < replab.gen.NiceIsomorphism
% Default implementation of a nice isomorphism for well-behaved group types

    properties (SetAccess = protected)
        sourceType % (`+replab.+gen.StaticFiniteGroupType`): Source type
    end

    methods

        function self = StaticNiceIsomorphism(sourceType, sourceGenerators, target)
        % Constructs a nice isomorphism for a static group type
        %
        % See `+replab.+gen.StaticFiniteGroupType.finishConstruction` for an explanation
        % of the overall mechanism.
        %
        % Args:
        %   sourceType (`.StaticFiniteGroupType`): Static group type for which we construct a nice isomorphism
        %   sourceGenerators (cell(1,\.*) of ``sourceType`` elements): Generators for this isomorphism source
        %   target (`+replab.FiniteGroup`) : Target of this isomorphism, with generators in one-to-one correspondance with the given ``sourceGenerators``
            self.source = sourceType.makeParentGroup(sourceGenerators, target, self);
            self.target = target;
            self.sourceType = sourceType;
            self.torusMap = [];
        end

        function setSource(self, source)
        % Sets the isomorphism source (internal)
        %
        % Hack used during construction of `+replab.AbstractGroup`, to solve a problem with
        % circular references.
        %
        % Args:
        %   source (`+replab.FiniteGroup`): Source to set
            self.source = source;
        end

    end

    methods % Implementations

        function l = sourceContains(self, s)
            l = true;
        end

        function t = imageElement(self, s)
            t = self.sourceType.imageElement(s);
        end

        function s = preimageElement(self, t)
            s = self.sourceType.preimageElement(t);
        end

    end

end
