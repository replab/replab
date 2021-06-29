classdef AbstractGroupRenamingIsomorphism < replab.FiniteIsomorphism

    methods

        function self = AbstractGroupRenamingIsomorphism(source, target)
            self.source = source;
            self.target = target;
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            t = self.target.imageLetters(self.source.factorizeLetters(s));
        end

        function s = preimageElement(self, t)
            s = self.source.imageLetters(self.target.factorizeLetters(t));
        end

    end

end
