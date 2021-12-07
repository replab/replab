classdef NormalCosets < replab.LeftCosets & replab.RightCosets
% Describes the set of normal cosets of a finite group
%
% Let $H$ be a normal subgroup of a group $G$.
% Then the left cosets are the sets $g H = \{ g h : h \in H \}$.
% Then the right cosets are the sets $H g = \{ h g : h \in H \}$.
% Because the subgroup is normal, we have $g H = H g$.

    methods % Implementation

        % Obj

        function l = laws(self)
            l = replab.laws.LeftCosetsLaws(self);
        end

        % Cosets

        function C = elements(self)
        % Returns the set of normal cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.NormalCoset`: Set of normal cosets
            C = cellfun(@(t) replab.NormalCoset(self.subgroup, t, self.group), self.transversal, 'uniform', 0);
        end

        function t = cosetRepresentative(self, g)
            t = cosetRepresentative@replab.LeftCosets(self, g);
        end

        function T = transversal(self)
            T = transversal@replab.RightCosets(self);
        end

    end

end
