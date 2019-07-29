classdef Monoid < replab.Domain
% Describes a monoid
    
    properties (SetAccess = protected)
        identity;
    end
    
    methods % Abstract methods
        
        function z = compose(self, x, y)
        % Returns the result of the binary operation applied to x and y
            f = self.composeFun;
            z = f(x, y);
        end

    end
    
    methods % Default implementations
        
        function z = composeAll(self, elements)
        % For self.composeAll({x1 x2 ... xn}), returns
        % x1*x2*...*xn (shown here in multiplaticative notation)
            if length(elements) == 0
                assert(isa(self, 'replab.Monoid'));
                z = self.identity;
            else
                z = elements{1};
                for i = 2:length(elements)
                    z = self.compose(z, elements{i});
                end
            end
        end
        
        function y = composeN(self, x, n)
        % Computes y = x^n by repeated squaring
        %
        % When "self" is a
        % - "Semigroup", we need n > 0
        % - "Monoid", we need n >= 0
        % - "Group", "n" is an arbitrary integer
            if n < 0
                assert(isa(self, 'replab.Group'));
                y = self.composeN(self.inverse(x), -n);
            elseif n == 0
                assert(isa(self, 'replab.Monoid'));
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
        
        function b = isIdentity(self, x)
        % Returns true if x is the identity, false otherwise
            b = self.eqv(x, self.identity);
        end
        
    end
end
