classdef IndexedFamilyLaws < replab.Laws
% Laws for indexed families

    properties (SetAccess = protected)
        I % `+replab.IndexedFamily`: Indexed family to test
        B % Domain of element indices as type vpi
        D % Domain of element indices as type double (capped by 2^53 -1)
    end

    methods

        function self = IndexedFamilyLaws(I)
        % Constructs a laws instance
        %
        % Args:
        %   I (`.IndexedFamily`): Indexed family to test
            self.I = I;
            self.B = replab.domain.vpi(1, I.size);
            if I.size > 2^53 - 1
                maxD = 2^53 - 1;
            else
                maxD = double(I.size);
            end
            self.D = replab.domain.intAsDouble(1, maxD);
        end

        function law_sample_is_contained_I(self, el)
        % Checks that every sample can be retrieved
            self.assert(~isempty(self.I.find(el)));
        end

        function law_at_find_roundtrip_B(self, ind)
        % Find/at roundtrip using the ``vpi`` index type
            a = self.I.at(ind);
            ind1 = self.I.find(a);
            self.B.assertEqv(ind, ind1);
        end

        function law_at_find_roundtrip_double_D(self, ind)
        % Find/at roundtrip using the ``double`` index type
            a = self.I.at(ind);
            ind1 = double(self.I.find(a));
            self.D.assertEqv(ind, ind1);
        end

        function law_at_find_roundtrip_string_B(self, ind)
        % Find/at roundtrip using the ``charstring`` index type
            str = strtrim(num2str(ind));
            a = self.I.at(str);
            ind1 = double(self.I.find(a));
            self.D.assertEqv(ind, ind1);
        end

    end

end
