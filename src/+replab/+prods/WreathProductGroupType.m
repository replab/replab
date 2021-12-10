classdef WreathProductGroupType < replab.gen.FiniteGroupType
% Describes the type of a wreath product of a permutation group acting on a finite group

    properties (SetAccess = protected)
        H % (`+replab.perm.PermutationGroupType`): Permutation group type
        A % (`+replab.FiniteGroupType`): Factor group type
        n % (integer): Number of copies of the factor group
    end

    methods

        function self = WreathProductGroupType(H, A)
            assert(isa(H, 'replab.perm.PermutationGroupType'));
            assert(isa(A, 'replab.FiniteGroupType'));
            self.H = H;
            self.A = A;
            n = H.domainSize;
            self.identity = {H.identity, repmat({A.identity}, 1, n)};
            self.n = n;
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = false; % guilty until proven innocent
            if ~self.H.eqv(x{1}, y{1})
                return
            end
            for i = 1:self.n
                if ~self.A.eqv(x{2}{i}, y{2}{i})
                    return
                end
            end
            b = true;
        end

        function g = sample(self)
            g = {self.H.sample, arrayfun(@(i) self.A.sample, 1:self.n, 'uniform', 0)};
        end

        % TotalOrder

        function c = compare(self, x, y)
            c = self.H.compare(x{1}, y{1});
            if c ~= 0
                return
            end
            for i = 1:self.n
                c = self.A.compare(x{2}{i}, y{2}{i});
                if c ~= 0
                    return
                end
            end
        end

        function z = compose(self, x, y)
            xh = x{1};
            xn = x{2};
            yh = y{1};
            yn = y{2};
            zh = self.H.compose(xh, yh);
            zn = arrayfun(@(i) self.A.compose(xn{yh(i)}, yn{i}), 1:self.n, 'uniform', 0);
            z = {zh, zn};
        end

        % Group

        function z = inverse(self, x)
            xh = x{1};
            xn = x{2};
            zh = self.H.inverse(xh);
            zn = cellfun(@(a) self.A.inverse(a), xn(zh), 'uniform', 0);
            z = {zh, zn};
        end

        % FiniteGroupTzpe

        function mu = niceIsomorphism(self, elements)
            nonIdentity = cell(1, 0);
            els = cellfun(@(el) el{2}, elements, 'uniform', 0);
            allA = horzcat(cell(1, 0), els{:});
            mask = cellfun(@(a) ~self.A.isIdentity(a), allA);
            A = self.A.groupWithGenerators(allA(mask));
            H = replab.S(self.n);
            mu = replab.prods.WreathProductNiceIsomorphism(self, H, A);
        end

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.prods.WreathProductGroupType') && self.H.isSameTypeAs(otherType.H) && self.A.isSameTypeAs(otherType.A);
        end

    end

end
