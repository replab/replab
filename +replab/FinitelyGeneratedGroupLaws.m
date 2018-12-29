classdef FinitelyGeneratedGroupLaws < replab.GroupLaws
    properties (SetAccess = protected)
        I; % index of generator
    end
    methods
        function self = FinitelyGeneratedGroupLaws(T)
            self = self@replab.GroupLaws(T);
            self.I = replab.domain.intAsDouble(1, T.nGenerators);
        end
    end
    methods
        function law_factorization_evaluateWord_T(self, t)
            w = self.T.factorization(t);
            t1 = self.T.evaluateWord(w);
            self.T.assertEqv(t, t1);
        end
        function law_generator_I(self, i)
            w = replab.Word.generator(i);
            t = self.T.generator(i);
            t1 = self.T.evaluateWord(w);
            self.T.assertEqv(t, t1);
        end
        function law_generatorInverse_I(self, i)
            t = self.T.inverse(self.T.generator(i));
            t1 = self.T.generatorInverse(i);
            self.T.assertEqv(t, t1);
        end
        function law_isTrivial(self)
            self.assert(self.T.isTrivial == (self.T.nGenerators == 0));
        end
        function law_freeGroup(self)
            self.assert(self.T.freeGroup.nGenerators == self.T.nGenerators);
        end
        function law_morphismFromFreeGroup_T(self, t)
            w = self.T.factorization(t);
            t1 = self.T.evaluateWord(w);
            f = self.T.morphismFromFreeGroup;
            t2 = f(w);
            self.T.assertEqv(t1, t2);
        end
        function morphismLaws = laws_morphismFromFreeGroup(self)
            morphismLaws = replab.GroupMorphismLaws(self.T.morphismFromFreeGroup, self.T.freeGroup, self.T);
        end
    end    
end
