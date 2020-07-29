classdef Isomorphism < replab.Morphism
% Describes an isomorphism between groups

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.IsomorphismLaws(self);
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

end
