classdef AbstractGroupNiceIsomorphism < replab.mrp.NiceIsomorphism

    methods

        function self = AbstractGroupNiceIsomorphism(source)
            self.source = source;
            self.target = source.permutationGroup;
        end

    end

    methods % Implementations

        function s = preimageElement(self, t)
            s = self.source.imageLetters(self.target.factorizeLetters(t));
        end


    end

end
