classdef RightCoset < replab.RightCoset

    methods

        function self = RightCoset(representative, subgroup, group)
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
            el = replab.bsgs.Cosets.rightRepresentative(self.subgroup.lexChain, el);
            b = self.type.eqv(self.representative, el);
        end

        function E = elementsSequence(self)
            H = self.subgroup.lexChain.allElements;
            g = self.representative;
            matrix = zeros(size(H));
            for i = 1:size(H, 2)
                matrix(:,i) = H(g,i);
            end
            matrix = sortrows(matrix')';
            E = replab.perm.Sequence(matrix);
        end

        function s = setProduct(self)
            s = replab.SetProduct(self.type, horzcat(self.subgroup.setProduct.sets, {{self.representative}}), false);
        end

    end

end
