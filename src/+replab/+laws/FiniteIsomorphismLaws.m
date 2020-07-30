classdef FiniteIsomorphismLaws < replab.laws.IsomorphismLaws & replab.laws.FiniteMorphismLaws

    methods

        function self = FiniteIsomorphismLaws(morphism)
            self@replab.laws.IsomorphismLaws(morphism);
            self@replab.laws.FiniteMorphismLaws(morphism);
        end

    end

end
