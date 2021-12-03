classdef NormalCosets < replab.perm.LeftCosets & replab.NormalCosets & replab.perm.RightCosets
% Describes the set of normal cosets of a finite group
%
% Let $H$ be a normal subgroup of a group $G$.
% Then the left cosets are the sets $g H = \{ g h : h \in H \}$.
% Then the right cosets are the sets $H g = \{ h g : h \in H \}$.
% Because the subgroup is normal, we have $g H = H g$.

    methods

        function self = NormalCosets(group, subgroup)
            self@replab.perm.LeftCosets(group, subgroup);
            self@replab.perm.RightCosets(group, subgroup);
        end

    end

    methods % Implementations

        function t = cosetRepresentative(self, g)
            t = cosetRepresentative@replab.perm.LeftCosets(self, g);
        end

        function T = transversal(self)
            T = transversal@replab.perm.RightCosets(self);
        end

    end

end
