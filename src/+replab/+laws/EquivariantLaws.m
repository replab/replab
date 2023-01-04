classdef EquivariantLaws < replab.laws.DomainLaws
% Law checks for a `+replab.Equivariant` space

    methods

        function self = EquivariantLaws(S)
            self@replab.laws.DomainLaws(S);
        end

        function law_project_twice_S(self, X)
        % Checks that an element is equivalent to itself
            [X1, err1] = self.S.project(X);
            [X2, err2] = self.S.project(X1);
            self.assertApproxEqual(X1, X2, err1+err2);
        end

    end

end
