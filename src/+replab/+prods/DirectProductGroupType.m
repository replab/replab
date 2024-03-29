classdef DirectProductGroupType < replab.gen.FiniteGroupType
% Describes the type of a direct product of group elements
%
% Note the duplication of the methods with `+replab.DirectProductGroup`, as we did not want to make the public-facing
% direct product group hierarchy too complex.

    properties (SetAccess = protected)
        factors % (cell(1,\*) of `+replab.FiniteGroupType`): Type of the factors
    end

    methods

        function self = DirectProductGroupType(factors)
            assert(all(cellfun(@(f) isa(f, 'replab.FiniteGroupType'), factors)));
            self.factors = factors;
            self.identity = cellfun(@(f) f.identity, factors, 'uniform', 0);
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = true;
            for i = 1:length(self.factors)
                if ~self.factors{i}.eqv(x{i}, y{i})
                    b = false;
                    return
                end
            end
        end

        function g = sample(self)
            g = cell(1, length(self.factors));
            for i = 1:length(self.factors)
                g{i} = self.factors{i}.sample;
            end
        end

        % TotalOrder

        function c = compare(self, x, y)
            for i = 1:length(self.factors)
                c = self.factors{i}.compare(x{i}, y{i});
                if c ~= 0
                    return
                end
            end
        end

        % Monoid

        function z = compose(self, x, y)
            z = cell(1, length(self.factors));
            for i = 1:length(self.factors)
                z{i} = self.factors{i}.compose(x{i}, y{i});
            end
        end

        % Group

        function xInv = inverse(self, x)
            xInv = cell(1, length(self.factors));
            for i = 1:length(self.factors)
                xInv{i} = self.factors{i}.inverse(x{i});
            end
        end

        % FiniteGroupType

        function mu = niceIsomorphism(self, elements)
            n = length(self.factors);
            factorGroups = cell(1, n);
            for i = 1:n
                factorElements = cell(1, 0);
                for j = 1:length(elements)
                    el = elements{j}{i};
                    if ~self.factors{i}.isIdentity(el)
                        factorElements{1,end+1} = el;
                    end
                end
                factorGroups{i} = self.factors{i}.groupWithGenerators(factorElements);
            end
            mu = replab.prods.DirectProductNiceIsomorphism(factorGroups);
        end

        function l = isSameTypeAs(self, otherType)
            l = false; % guilty until proven innocent
            if ~isa(otherType, 'replab.prods.DirectProductGroupType')
                return
            end
            if length(self.factors) ~= length(otherType.factors)
                return
            end
            for i = 1:length(self.factors)
                if ~self.factors{i}.isSameTypeAs(otherType.factors{i})
                    return
                end
            end
            l = true;
        end

    end

end
