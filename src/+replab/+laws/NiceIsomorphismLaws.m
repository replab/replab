classdef NiceIsomorphismLaws < replab.laws.FiniteIsomorphismLaws

    methods

        function self = NiceIsomorphismLaws(morphism)
            self@replab.laws.FiniteIsomorphismLaws(morphism);
        end

    end

    methods % Laws

        function law_source_must_be_generic_(self)
            self.assert(isa(self.morphism.source.type, 'replab.gen.FiniteGroupType'));
        end

    end

end
