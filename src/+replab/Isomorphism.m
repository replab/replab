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

        function m = lambda(source, target, preimageElementFun, imageElementFun)
        % Creates a morphism from an image function
        %
        % Args:
        %   source (`.Group`): Source group
        %   target (`.Group`): Target group
        %   preimageElementFun (function_handle): Function computing preimages of elements
        %   imageElementFun (function_handle): Function computing images of elements
        %
        % Returns:
        %   `+replab.Isoorphism`: Constructed isomorphism
            m = replab.mrp.LambdaIsomorphism(source, target, preimageElementFun, imageElementFun);
        end

    end

end
