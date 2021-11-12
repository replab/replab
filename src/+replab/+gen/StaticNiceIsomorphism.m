classdef StaticNiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        sourceType % (`+replab.+gen.StaticFiniteGroupType`): Source type
    end

    methods

        function self = StaticNiceIsomorphism(sourceType, sourceGenerators, target)
        % Constructs a nice isomorphism for a static group type
        %
        % Args:
        %   type (`.StaticFiniteGroupType`): Static group type
            self.source = sourceType.makeIsomorphismSourceGroup(sourceGenerators, target, self);
            self.target = target;
            self.sourceType = sourceType;
            self.torusMap = [];
        end

        function setSource(self, source)
        % Sets the isomorphism source (internal)
        %
        % Hack used during construction of `+replab.AbstractGroup`, as the group is constructed
        % before the group type.
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
