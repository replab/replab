classdef Monoid < replab.cat.Domain
% Defines a group
    
    properties (SetAccess = protected)
        identity;
    end
    
    methods % Abstract methods
        
        function z = compose(self, x, y)
        % Returns the result of the group binary operation applied
        % to x and y
            f = self.composeFun;
            z = f(x, y);
        end
        
    end

    methods % Methods with default implementations
        
        function b = isIdentity(self, x)
        % Returns true if x is the identity, false otherwise
            b = self.eqv(x, self.identity);
        end
        
        function z = composeMany(self, varargin)
        % For self.composeMany(x1, x2, ... xn), returns
        % x1*x2*...*xn (shown here in multiplaticative notation)
            if length(varargin) == 0
                z = self.identity;
            else
                z = varargin{1};
                for i = 2:length(varargin)
                    z = self.compose(z, varargin{i});
                end
            end
        end
        
        function y = composeN(self, x, n)
        % Computes y = x^n by repeated squaring
            if n == 0
                y = self.identity;
            elseif n == 1
                y = x;
            elseif n == 2
                y = self.compose(x, x);
            else
                y = self.identity;
                while n > 1
                    if mod(n, 2) == 0 % n even
                        n = n / 2;
                    else % n odd
                        y = self.compose(x, y);
                        n = (n - 1)/2;
                    end
                    x = self.compose(x, x);
                end
                y = self.compose(x, y);
            end
        end
        
    end
    
    methods % Laws
        
        function law_associativity_DDD(self, x, y, z)
        % Checks associativity of group binary operation
            xy = self.compose(x, y);
            yz = self.compose(y, z);
            self.assertEqv(self.compose(xy, z), self.compose(x, yz));
        end

        function law_composeN_DN010(self, x, n)
            xn1 = self.identity;
            for i = 1:n
                xn1 = self.compose(xn1, x);
            end
            if n > 0
                xn2 = self.compose(self.composeN(x, n - 1), x);
                self.assertEqv(xn1, xn2);
            end
            xn3 = self.composeN(x, n);
            self.assertEqv(xn1, xn3);
        end

        function law_identity(self)
            self.assertTrue(self.isIdentity(self.identity));
        end

    end

end
