classdef NormalCoset < replab.NormalCoset & replab.PermutationFiniteSet

    methods

        function self = NormalCoset(representative, subgroup, group)
            assert(group.hasSameTypeAs(subgroup));
            self.type = group.type;
            self.representative_ = representative;
            self.group = group;
            self.subgroup = subgroup;
        end

    end

    methods % Implementations

        function b = contains(self, el)
            el = replab.bsgs.Cosets.leftRepresentative(self.subgroup.lexChain, el); % same implementation as `.LeftCoset`
            b = self.group.eqv(self.representative, el);
        end

        function E = elementsSequence(self)
            H = self.subgroup.lexChain.allElements; % same implementation as `.LeftCoset`
            g = self.representative;
            matrix = sortrows(g(H'))';
            E = replab.perm.Sequence(matrix);
        end

    end

end
