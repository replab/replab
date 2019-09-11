classdef Group < replab.Group
    
    properties (SetAccess = protected)
        header;
        eqvFun;
        sampleFun;
        composeFun;        
        inverseFun;
    end
    
    methods
        
        function self = Group(header, eqvFun, sampleFun, ... % Domain
                              composeFun, identity, ... % Monoid
                              inverseFun) % Group
            self.header = header;
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
            self.composeFun = composeFun;
            self.identity = identity;
            self.inverseFun = inverseFun;
        end
        
        function str = headerStr(self)
            str = self.header;
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Group(self), ...
                {'header'} ...
                );
        end
        
        % Domain methods
        
        function b = eqv(self, t, u)
            f = self.eqvFun;
            b = f(t, u);
        end

        function t = sample(self)
            f = self.sampleFun;
            t = f();
        end

        % Monoid methods
        
        function z = compose(self, x, y)
            f = self.composeFun;
            z = f(x, y);
        end
        
        % Group methods
        
        function xInv = inverse(self, x)
            f = self.inverseFun;
            xInv = f(x);
        end
        
    end
    
end
