classdef RightCoset < replab.Coset
% Describes a right coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The right coset is $H g = \{ h g : h \in H\}$.

    methods % Implementations

        % Domain

        function s = sample(self)
            s = self.type.compose(self.subgroup.sample, self.representative);
        end

        % FiniteSet

        function s = setProduct(self)
            s = replab.SetProduct(self.type, horzcat(self.subgroup.setProduct.sets, {{self.representative}}), false);
        end

    end

end
