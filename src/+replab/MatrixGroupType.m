classdef MatrixGroupType < replab.gen.FiniteGroupType

    properties (SetAccess = protected)
        n % (integer): Matrix size
    end

    methods

        function self = MatrixGroupType(n)
            self.n = n;
            self.identity = replab.cyclotomic.eye(n);
        end

    end

    methods % Implementations

        function G = groupWithGenerators(self, generators, varargin)
            error('Not implemented');
        end

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.MatrixGroupType') && self.n == otherType.n;
        end

        % Domain

        function l = eqv(self, x, y)
            l = isequal(x, y);
        end

        % TotalOrder

        function c = compare(self, x, y)
            xisid = all(all(x == eye(self.n)));
            yisid = all(all(y == eye(self.n)));
            if xisid && yisid
                c = 0;
            elseif xisid && ~yisid
                c = -1;
            elseif ~xisid && yisid
                c = 1;
            else
                x = double(x).';
                y = double(y).';
                x = x(:);
                y = y(:);
                x = [real(x).'; imag(x).'];
                y = [real(y).'; imag(y).'];
                v = x - y;
                ind = find(v ~= 0, 1);
                c = [sign(v(ind)) 0];
                c = c(1);
            end
        end

        % Group

        function z = compose(self, x, y)
            z = x * y;
        end

        function xInv = inverse(self, x)
            xInv = inv(x);
        end


    end

end
