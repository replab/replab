classdef FiniteIsomorphismLaws < replab.IsomorphismLaws & replab.FiniteMorphismLaws

    methods

        function self = FiniteIsomorphismLaws(morphism)
            self@replab.IsomorphismLaws(morphism);
            self@replab.FiniteMorphismLaws(morphism);
        end

    end

end
