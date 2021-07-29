classdef EquivariantLaws < replab.laws.DomainLaws
% Law checks for a `+replab.Equivariant` space

    methods

        function self = EquivariantLaws(T)
            self@replab.laws.DomainLaws(T);
        end

        function law_project_twice_T(self, X)
        % Checks that an element is equivalent to itself
            [X1, err] = self.T.project(X);
            self.assertApproxEqual(X, X1, err);
        end

    end

end
