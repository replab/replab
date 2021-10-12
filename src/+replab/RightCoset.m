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

        % Coset

        function [l, r] = factorizeShortRepresentativeLetters(self)
        % Returns a tentatively short word corresponding to an element of this coset
        %
        % An effort is made to identify a short word, but without optimality guarantees.
        %
        % Returns
        % -------
        %   l: integer(1,\*)
        %     Letters of the word representing an element of ``self``
        %   r: element of `.group`
        %     Represented coset element
            [linv, rinv] = self.subgroup.leftCoset(self.type.inverse(self.representative), 'group', self.group).factorizeShortRepresentativeLetters;
            l = replab.fp.Letters.inverse(linv);
            r = self.type.inverse(rinv);
        end

    end

end
