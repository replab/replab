classdef NormalCosets < replab.LeftCosets & replab.RightCosets
% Describes the set of normal cosets of a finite group
%
% Let $H$ be a normal subgroup of a group $G$.
% Then the left cosets are the sets $g H = \{ g h : h \in H \}$.
% Then the right cosets are the sets $H g = \{ h g : h \in H \}$.
% Because the subgroup is normal, we have $g H = H g$.

    methods

        function self = NormalCosets(group, subgroup)
            self@replab.LeftCosets(group, subgroup);
            self@replab.RightCosets(group, subgroup);
        end

        function s = cardinality(self)
            s = size@replab.LeftCosets(self);
        end

        function C = coset(self, g)
        % Returns the normal coset containing the given element
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   `+replab.NormalCoset`: Normal coset
            C = replab.NormalCoset(self.group, self.subgroup, self.cosetRepresentative(g));
        end

        function t = cosetRepresentative(self, g)
            t = cosetRepresentative@replab.LeftCosets(self, g);
        end

        function T = transversal(self)
            T = transversal@replab.RightCosets(self);
        end

        function C = elements(self)
        % Returns the set of normal cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.NormalCoset`: Set of normal cosets
            C = cellfun(@(t) replab.NormalCoset(self.group, self.subgroup, t), self.transversal, 'uniform', 0);
        end

    end

end
