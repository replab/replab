classdef EquivariantLaws < replab.laws.DomainLaws
% Law checks for a `+replab.Equivariant` space

    methods

        function self = EquivariantLaws(S)
            self@replab.laws.DomainLaws(S);
        end

        function law_project_twice_S(self, X)
        % Checks that an element is equivalent to itself
            [X1, err] = self.S.project(X);
            self.assertApproxEqual(X, X1, err);
        end

    end

end
