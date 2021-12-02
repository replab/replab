classdef LeftCoset < replab.LeftCoset

    methods

        function self = LeftCoset(representative, subgroup, group)
            assert(group.hasSameTypeAs(subgroup));
            self.type = group.type;
            self.representative_ = representative;
            self.group = group;
            self.subgroup = subgroup;
        end

    end

    methods % Implementations

        % FiniteSet

        function b = contains(self, el)
            el = replab.bsgs.Cosets.leftRepresentative(self.subgroup.lexChain, el);
            b = self.type.eqv(self.representative, el);
        end

        function E = elementsSequence(self)
            H = self.subgroup.lexChain.allElements;
            g = self.representative;
            matrix = sortrows(g(H'))';
            E = replab.perm.Sequence(matrix);
        end

        function s = setProduct(self)
            s = replab.SetProduct(self.type, horzcat({{self.representative}}, self.subgroup.setProduct.sets), false);
        end

    end

end
