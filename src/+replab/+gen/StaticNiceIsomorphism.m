classdef StaticNiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        sourceType % (`+replab.+gen.StaticFiniteGroupType`): Source type
    end

    methods

        function self = StaticNiceIsomorphism(type)
        % Constructs a nice isomorphism for a static group type
        %
        % Args:
        %   type (`.StaticFiniteGroupType`): Static group type
            self.sourceType = type;
            sourceGenerators = type.sourceGenerators;
            sourceFun = @(iso) type.makeSource(sourceGenerators, iso);
            self.finishConstruction(sourceFun, sourceGenerators, type.niceType);
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
