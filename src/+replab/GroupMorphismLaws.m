classdef GroupMorphismLaws < replab.Laws
    
    properties (SetAccess = protected)
        fun; % morphism function handle
        S; % source group
        T; % target group
    end
    
    methods
        
        function self = GroupMorphismLaws(fun, S, T)
            self.fun = fun;
            self.S = S;
            self.T = T;
        end
        
    end
    
    methods % LAWS
        
        function law_inverse_S(self, s)
            t = self.fun(s);
            sI = self.S.inverse(s);
            tI1 = self.T.inverse(t);
            tI2 = self.fun(sI);
            self.T.assertEqv(tI1, tI2);
        end
        
        function law_composition_SS(self, s1, s2)
            s12 = self.S.compose(s1, s2);
            t1 = self.fun(s1);
            t2 = self.fun(s2);
            t12_1 = self.fun(s12);
            t12_2 = self.T.compose(t1, t2);
            self.T.assertEqv(t12_1, t12_2);
        end
        
        function law_identity(self)
            self.T.assertEqv(self.T.identity, self.fun(self.S.identity));
        end
        
    end
    
end
