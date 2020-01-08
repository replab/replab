classdef Monoid < replab.Domain
% Describes a monoid
%
% See https://en.wikipedia.org/wiki/Monoid

    properties (SetAccess = protected)
        identity % Monoid identity element
    end

    methods % Abstract methods

        function z = compose(self, x, y)
        % Composes two monoid/group elements
        %
        % Args:
        %   x (element): Left hand side of the binary operation
        %   y (element): Right hand side of the binary operation
        %
        % Returns:
        %   element: Result of the binary operation ``x op y``
            error('Abstract');
        end

    end

    methods % Default implementations

        function z = composeAll(self, elements)
        % Composes a sequence of monoid elements
        %
        % For ``self.composeAll({x1 x2 ... xn})``, returns
        % ``x1*x2*...*xn`` (shown here in multiplaticative notation)
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
        % Computes ``y = x^n`` by repeated squaring
        %
        % When "self" is a
        % - `.Semigroup`, we need n > 0
        % - `.Monoid`, we need n >= 0
        % - `.Group`, ``n`` is an arbitrary integer
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

    methods (Static)

        function monoid = lambda(header, eqvFun, sampleFun, composeFun, identity)
        % Constructs a monoid from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the `eqv` method
        %   sampleFun (function_handle): Handle implementing the `sample` method
        %   composeFun (function_handle): Handle implementing the `compose` method
        %   identity (element): Identity element of this monoid
        %
        % Returns:
        %   replab.Monoid: The constructed monoid
            monoid = replab.lambda.Monoid(header, eqvFun, sampleFun, composeFun, identity);
        end

    end

end
