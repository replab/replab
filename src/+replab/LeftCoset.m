classdef LeftCoset < replab.Coset
% Describes a left coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The left coset is $g H = \{ g h : h \in H\}$.

    methods % Implementations

        % Domain

        function s = sample(self)
            s = self.type.compose(self.representative, self.subgroup.sample);
        end

        % FiniteSet

        function s = setProduct(self)
            s = replab.SetProduct(self.type, horzcat({{self.representative}}, self.subgroup.setProduct.sets), false);
        end

    end

end
