classdef ConjugacyClass < replab.FiniteSet
% Describes a conjugacy class of a finite group
%
% A conjugacy class containing the representative $r \in G$ is the set $\{g r g^{-1} : g \in G \}$.
%
% The centralizer of $r$ in $G$ is the subgroup $C_{G}(r) = \{ g r g^{-1} == r : g \in G \}$.
%
% Thus, the left cosets $G/C_{G}(r) = \{ g C_{G}(r) : g \in G \}$ are in one to one correspondence with
% the elements of the conjugacy class.

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group containing this conjugacy class
    end

    methods

        function l = knownRepresentativeCentralizer(self)
        % Returns whether the centralizer of `.representative` is known
        %
        % Returns:
        %   logical: True if the centralizer has already been computed
            error('Abstract');
        end

        function G = representativeCentralizer(self)
        % Returns the centralizer of `.representative` in `.group`
        %
        % Returns:
        %   `.FiniteGroup`: Representative centralizer
            error('Abstract');
        end

    end

    methods

        function o = elementOrder(self)
        % Returns the order of the elements in this conjugacy class
        %
        % Returns:
        %   integer: Element order
            error('Abstract');
        end

        function c1 = imap(self, f)
        % Maps this conjugacy class under an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %
        % Returns:
        %   `.ConjugacyClass`: The conjugacy class mapped under ``f``, expressed as a subset of ``f.image``
            group1 = self.group.imap(f);
            rep1 = f.imageElement(self.representative);
            args = {};
            if self.knownRepresentativeCentralizer
                rc1 = self.representativeCentralizer.imap(f);
                args = {'centralizer', rc1};
            end
            c1 = group1.conjugacyClass(rep1, 'isCanonical', f.preservesTypeOrder, args{:});
        end

    end

    methods % Implementations

        % Str

        function s = shortStr(self, maxColumns)
            s = sprintf('ConjugacyClass of %s in %s', replab.shortStr(self.representative, maxColumns), replab.shortStr(self.group, maxColumns));
            if length(s) > maxColumns
                s = sprintf('ConjugacyClass of %s', replab.shortStr(self.representative, maxColumns));
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.ConjugacyClassLaws(self);
        end

        % Domain

        function s = sample(self)
            t = self.group.sample;
            s = self.group.leftConjugate(t, self.representative);
        end

    end

end
