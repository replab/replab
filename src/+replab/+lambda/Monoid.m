classdef Monoid < replab.Monoid
    
    properties (SetAccess = protected)
        header;
        eqvFun;
        sampleFun;
        composeFun;
    end

    methods
        
        function self = Monoid(header, eqvFun, sampleFun, ... % Domain
                               composeFun, identity) % Monoid
            self.header = header;
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
            self.composeFun = composeFun;
            self.identity = identity;
        end
        
        function str = headerStr(self)
            str = self.header;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Monoid(self), ...
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
            error('Not implemented');
            f = self.composeFun;
            z = f(x, y);
        end
        
    end
    
end
