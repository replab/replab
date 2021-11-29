classdef TotalOrderLaws < replab.laws.DomainLaws

    methods

        function self = TotalOrderLaws(S)
            self@replab.laws.DomainLaws(S);
        end

    end

    methods % Laws

        % compare

        function law_reflexive_S(self, x)
            self.assert(self.S.compare(x, x) == 0);
        end

        function law_transitive_SSS(self, x, y, z)
            if self.S.compare(x, y) <= 0 && self.S.compare(y, z) <= 0
                self.assert(self.S.compare(x, z) <= 0);
            end
        end

        % sort

        function law_sort_SSSS(self, s1, s2, s3, s4)
            array = [s1 s2 s3 s4];
            I = self.S.sort(array);
            array = array(I);
            for i = 1:4
                for j = i:4
                    self.assert(self.S.compare(array(i), array(j)) <= 0);
                end
            end
        end

    end

end
