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

        function s = nElements(self)
            s = nElements@replab.LeftCosets(self);
        end

        function t = cosetRepresentative(self, g)
            t = cosetRepresentative@replab.LeftCosets(self, g);
        end

        function T = transversal(self)
            T = transversal@replab.RightCosets(self);
        end

        function T = computeTransversal(self)
            T = computeTransversal@replab.RightCosets(self);
        end

        function C = elements(self)
        % Returns the set of normal cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.NormalCoset`: Set of normal cosets
            C = cellfun(@(t) replab.NormalCoset(self.subgroup, t, self.group), self.transversal, 'uniform', 0);
        end

    end

end
