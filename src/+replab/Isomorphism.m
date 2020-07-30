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
            I = replab.fm.Inverse(self);
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
            m = replab.fm.Identity(group);
        end

    end

end
