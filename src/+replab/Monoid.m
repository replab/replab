classdef Monoid < replab.Domain
% Describes a monoid
%
% See `Monoid on Wikipedia <https://en.wikipedia.org/wiki/Monoid>`_

    properties (SetAccess = protected)
        identity % (element): Monoid identity element
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

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.MonoidLaws(self);
        end

    end

    methods % Monoid operations

        function z = composeAll(self, elements)
        % Composes a sequence of monoid elements
        %
        % For ``self.composeAll({x1 x2 ... xn})``, returns ``x1*x2*...*xn``
        % (shown here in multiplaticative notation)
        %
        % Args:
        %   elements (cell(1,\*) of element): Elements to compose
        %
        % Returns:
        %   element: Composition of the given elements
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
        % Computes the power of a given element by repeated squaring
        %
        % Args:
        %   x (element): Monoid/group element
        %   n (integer): Power; if this is not a group, we need ``n >= 0``
        %
        % Returns:
        %   element: The element ``x^n``
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
        % Tests if the given element is the identity
        %
        % Args:
        %   x (element): Monoid/Group element
        %
        % Returns:
        %   logical: True if ``x`` is the identity, false otherwise
            b = self.eqv(x, self.identity);
        end

        function b = areCommuting(self, x, y)
        % Returns whether the two given elements commute
        %
        % Args:
        %   x (element): Monoid/Group element
        %   y (element): Monoid/Group element
        %
        % Returns:
        %   logical: True if ``x y = y x``, false otherwise
            b = self.eqv(self.compose(x, y), self.compose(y, x));
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
        %   `.Monoid`: The constructed monoid
            monoid = replab.lambda.Monoid(header, eqvFun, sampleFun, composeFun, identity);
        end

    end

end
