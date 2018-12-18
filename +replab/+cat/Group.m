classdef Group < replab.cat.Domain
% Defines a group
    
    properties (SetAccess = protected)
        identity;
    end
    
    methods % Abstract methods
        
        %ABSTRACT z = compose(self, x, y)
        %ABSTRACT xInv = inverse(self, x)
        
    end

    methods % Methods with default implementations
        
        function b = isIdentity(self, x)
            b = self.eqv(x, self.identity);
        end
        
        function x = conjugate(self, by, on)
        % Returns "on" conjugated by "by", i.e.
        % by * on * by.inverse
            x = self.composeWithInverse(self.compose(by, on), by);
        end
        
        function z = composeWithInverse(self, x, y)
            z = self.compose(x, self.inverse(y));
        end
        
        function z = composeMany(self, varargin)
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
        % Computes y = x^e by repeated squaring
            if n < 0
                y = self.composeN(self.inverse(x), -n);
            elseif n == 0
                y = self.identity;
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

        function law_composeN_DZ10(self, x, n)
            xn1 = self.identity;
            if n < 0
                for i = 1:(-n)
                    xn1 = self.composeWithInverse(xn1, x);
                end
            else
                for i = 1:n
                    xn1 = self.compose(xn1, x);
                end
            end
            xn2 = self.compose(self.composeN(x, n - 1), x);
            xn3 = self.composeN(x, n);
            self.assertEqv(xn1, xn2);
            self.assertEqv(xn1, xn3);
        end

        function law_identity(self)
            self.assertTrue(self.isIdentity(self.identity));
        end

        function law_inverse_D(self, x)
            xI = self.inverse(x);
            id1 = self.compose(x, xI);
            id2 = self.compose(xI, x);
            self.assertTrue(self.isIdentity(id1));
            self.assertTrue(self.isIdentity(id2));
        end
        
        function law_inverse_compatible_with_compose_DD(self, x, y)
            xy = self.compose(x, y);
            yIxI = self.compose(self.inverse(y), self.inverse(x));
            self.assertEqv(self.inverse(xy), yIxI);
        end
                
    end

end
