classdef SequenceLaws < replab.Laws
% Laws for indexed families

    properties (SetAccess = protected)
        S % (`+replab.Sequence`): Sequence to test
        B % Domain of element indices as type vpi
        D % Domain of element indices as type double (capped by 2^53 -1)
    end

    methods

        function self = SequenceLaws(S)
        % Constructs a laws instance
        %
        % Args:
        %   I (`+replab.Sequence`): Indexed family to test
            self.S = S;
            self.B = replab.domain.vpi(1, S.nElements);
            if S.nElements > 2^53 - 1
                maxD = 2^53 - 1;
            else
                maxD = double(S.nElements);
            end
            self.D = replab.domain.intAsDouble(1, maxD);
        end

        function law_sample_is_contained_S(self, el)
        % Checks that every sample can be retrieved
            self.assert(~isempty(self.S.find(el)));
        end

        function law_at_find_roundtrip_B(self, ind)
        % Find/at roundtrip using the ``vpi`` index type
            a = self.S.at(ind);
            ind1 = self.S.find(a);
            self.B.assertEqv(ind, ind1);
        end

        function law_at_find_roundtrip_double_D(self, ind)
        % Find/at roundtrip using the ``double`` index type
            a = self.S.at(ind);
            ind1 = double(self.S.find(a));
            self.D.assertEqv(ind, ind1);
        end

        function law_at_find_roundtrip_string_B(self, ind)
        % Find/at roundtrip using the ``charstring`` index type
            str = strtrim(num2str(ind));
            a = self.S.at(str);
            ind1 = double(self.S.find(a));
            self.D.assertEqv(ind, ind1);
        end

    end

end
