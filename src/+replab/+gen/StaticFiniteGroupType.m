classdef StaticFiniteGroupType < replab.gen.FiniteGroupType

    properties (SetAccess = protected)
        niceType % (`+replab.FiniteGroupType`): Target type
    end

    methods % Internal methods

        function t = imageElement(self, s)
        % Computes the (nice) image of an element of this type
        %
        % Args:
        %   s (element of ``self``): Element of this type
        %
        % Returns:
        %   element of `.targetType`: Image of the element
            error('Abstract');
        end

        function mu = isomorphism(self)
            mu = self.cached('isomorphism', @() replab.gen.StaticNiceIsomorphism(self));
        end

        function s = preimageElement(self, t)
        % Computes the preimage of an element
        %
        % Args:
        %   t (element of `.isomorphism` target): Nice image of an element
        %
        % Returns:
        %   element of ``self``: Preimage of the element
            error('Abstract');
        end

        function S = makeSource(self, generators, niceIsomorphism)
        % Constructs the group that contains all type elements
            S = replab.gen.FiniteGroup(self, generators, 'niceIsomorphism', niceIsomorphism);
        end

    end

    methods

        function iso = constructIsomorphism(self, elements)
            iso = self.isomorphism;
        end

    end

end
