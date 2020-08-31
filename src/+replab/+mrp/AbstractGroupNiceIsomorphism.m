classdef AbstractGroupNiceIsomorphism < replab.FiniteIsomorphism

    methods

        function self = AbstractGroupNiceIsomorphism(source)
            self.source = source;
            self.target = source.permutationGroup;
        end

    end

    methods % Implementations


        function T = imageGroup(self, S)
            T = S.niceGroup;
        end

        function t = imageElement(self, s)
            t = self.source.niceImage(s);
        end

        function s = preimageElement(self, t)
            s = self.source.imageLetters(self.target.factorizeLetters(t));
        end

        function S = preimageGroup(self, T)
            gens = cellfun(@(t) self.preimageElement(t), T.generators, 'uniform', 0);
            S = self.source.niceSubgroup(gens, T.order, T);
        end


    end

end
