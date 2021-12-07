classdef Isomorphism < replab.Morphism
% Describes an isomorphism between groups

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.IsomorphismLaws(self);
        end

    end

    methods (Access = protected)

        function I = computeInverse(self)
            I = replab.mrp.Inverse(self);
        end

    end


    methods % Morphism composition

        function I = inverse(self)
            I = self.cached('inverse', @() self.computeInverse);
        end

    end

    methods % Preimages

        function s = preimageElement(t)
        % Returns an the preimage of an element of the range
        %
        % Returns the ``s`` such that ``self.imageElement(s) == t`` .
        %
        % If the argument is not in the range of this isomorphism, the behavior is undefined.
        %
        % Args:
        %   t (element of `.target`): Element to compute the preimage of
        %
        % Returns:
        %   element of `.source`: Preimage
            error('Abstract');
        end

    end

    methods (Static)

        function m = identity(group)
        % Returns the identity morphism from a group to itself
        %
        % Args:
        %   group (`.Group`): Group
        %
        % Returns:
        %   `.Isomorphism`: The identity automorphism on the given group
            m = replab.mrp.Identity(group);
        end

        function m = lambda(source, target, preimageElementFun, imageElementFun, torusMap)
        % Creates a morphism from an image function
        %
        % Args:
        %   source (`.Group`): Source group
        %   target (`.Group`): Target group
        %   preimageElementFun (function_handle): Function computing preimages of elements
        %   imageElementFun (function_handle): Function computing images of elements
        %   torusMap (integer(\*,\*)): Torus map to use in the morphism construction
        %
        % Returns:
        %   `+replab.Isomorphism`: Constructed isomorphism
            if nargin < 5
                torusMap = [];
            end
            m = replab.mrp.LambdaIsomorphism(source, target, preimageElementFun, imageElementFun, torusMap);
        end

    end

end
