classdef DoubleCosets < replab.Obj
% Describes the set of double cosets in a group
%
% A double coset is a set of the form ``{ H g K } = { h g k : h \in H, k \in K}`` for subgroups
% ``H`` and ``K`` of a group ``G``

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group
        leftSubgroup % (`.FiniteGroup`): Subgroup of `.group`
        rightSubgroup % (`.FiniteGroup`): Subgroup of `.group`
    end

    methods

        function t = cosetRepresentative(self, g)
        % Returns the canonical coset representative corresponding to the given element
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   t (element of `.group`): Double coset canonical representative
            error('Abstract');
        end

        function s = nElements(self)
        % Returns the number of double cosets
        %
        % Returns:
        %   integer: Number of double cosets
            s = length(self.transversal);
        end

        function C = elements(self)
        % Returns the set of double cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.DoubleCoset`: Set of double cosets
            C = cellfun(@(t) self.leftSubgroup.doubleCoset(t, self.rightSubgroup, 'group', self.group), self.transversal, 'uniform', 0);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            error('Abstract');
        end

    end

end
