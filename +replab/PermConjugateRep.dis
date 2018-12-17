classdef PermConjugateRep < replab.Rep
    
    properties (SetAccess = private)
        conjugateBy;
        child;
    end
    
    methods
        
        function self = PermConjugateRep(conjugateBy, child)
            self.group = child.group;
            self.d = child.d;
            self.field = child.field;
            self.conjugateBy = conjugateBy;
            self.child = child;
        end
        
        function prettyDisp(self, spaces)
            disp(sprintf('%s Representation conjugated by %s:', spaces, num2str(self.conjugateBy)));
            self.child.prettyDisp([spaces '  ']);
        end
                
        function M = image(self, permutation)
            M = self.child.image(permutation);
            M(self.conjugateBy, self.conjugateBy) = M;
        end
        
    end
        
end
