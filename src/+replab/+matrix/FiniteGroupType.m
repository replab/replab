classdef FiniteGroupType < replab.gen.FiniteGroupType

    properties (SetAccess = protected)
        d % (integer): Matrix size
    end

    methods

        function self = FiniteGroupType(d)
            self.d = d;
            self.identity = replab.cyclotomic.eye(d);
        end

    end

    methods % Implementations

        % Domain

        function l = eqv(self, x, y)
            l = isequal(x, y);
        end

        % TotalOrder

        function c = compare(self, x, y)
            e = self.identity;
            x_diff = sum(sum(x ~= e));
            y_diff = sum(sum(y ~= e));
            if x_diff < y_diff
                c = -1;
                return
            elseif x_diff > y_diff
                c = 1;
                return
            end
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

        % Group

        function z = compose(self, x, y)
            z = x * y;
        end

        function xInv = inverse(self, x)
            xInv = inv(x);
        end

        % FiniteGroupType

        function mu = niceIsomorphism(self, elements)
            mask = cellfun(@(g) self.isIdentity(g), elements);
            mu = replab.matrix.NiceIsomorphism(self.d, elements(~mask));
        end

        function G = groupWithGenerators(self, generators, varargin)
            generators = cellfun(@(g) replab.cyclotomic(g), generators, 'uniform', 0);
            G = replab.MatrixGroup(self.d, generators, 'type', self, varargin{:});
        end

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.matrix.FiniteGroupType') && self.d == otherType.d;
        end

        function G = subgroupWithGenerators(self, group, generators, varargin)
            generators = cellfun(@(g) replab.cyclotomic(g), generators, 'uniform', 0);
            G = replab.MatrixGroup(self.d, generators, 'type', self, 'niceIsomorphism', group.niceIsomorphism, varargin{:});
        end

    end

end
