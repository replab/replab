classdef Group < replab.Monoid
% Describes a group
    
    methods % Abstract methods
        
        function xInv = inverse(self, x)
        % Returns the inverse element xInv of x, such that
        % x xInv = identity
            f = self.inverseFun;
            xInv = f(x);
        end
        
    end

    methods % Methods with default implementations
        
        function x = leftConjugate(self, by, on)
        % Returns "on" conjugated by "by", that is
        % x = by * on * by^-1 in multiplicative notation
            x = self.composeWithInverse(self.compose(by, on), by);
        end
        
        function z = composeWithInverse(self, x, y)
        % Returns x * y^-1 (shown here in multiplicative notation)
            z = self.compose(x, self.inverse(y));
        end
        
    end

    methods (Static)
        
        function group = lambda(header, eqvFun, sampleFun, ...
                                composeFun, identity, inverseFun)
            group = replab.lambda.Group(header, eqvFun, sampleFun, ...
                                        composeFun, identity, inverseFun);
        end
        
    end
    
end
