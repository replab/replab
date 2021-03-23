classdef Coset < replab.FiniteSet

    properties (SetAccess = protected)
        parent % (`.FiniteGroup`): Group containing this coset
        group % (`.FiniteGroup`): Group
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism to a permutation group
        groupChain % (`+replab.+bsgs.Chain`): Group chain with base in lexicographic order
    end

    methods

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
            error('Abstract');
        end

        function [w, r] = factorizeShortRepresentativeWord(self)
        % Returns a tentatively short word corresponding to an element of this coset
        %
        % An effort is made to identify a short word, but without optimality guarantees.
        %
        % Returns
        % -------
        %   w: charstring
        %     Word representing an element of ``self``
        %   r: element of `.group`
        %     Represented coset element
            [l, r] = self.factorizeShortRepresentativeLetters;
            w = replab.fp.Letters.print(l, self.group.names);
        end

    end

    methods % Implementations

        % Domain

        function l = laws(self)
            l = replab.laws.CosetLaws(self);
        end

        function b = eqv(self, lhs, rhs)
            b = self.parent.eqv(lhs, rhs);
        end

    end

end
