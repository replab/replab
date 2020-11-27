classdef TensorRep < replab.Rep
% A tensor product of representations
%
% All factor representations must be defined on the same group

    properties
        factors % (cell(1,\*) of `+replab.Rep`): Factor representations
    end

    methods

        function self = TensorRep(group, field, factors)
        % Constructs a tensor representation from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   group (`+replab.CompactGroup`): Common group
        %   field ({'R', 'C'}): Real or complex field
        %   blocks (cell(1,\*) of `+replab.Rep`): Factor representations
            replab.rep.assertCompatibleFactors(group, field, factors);
            factorsAllUnitary = cellfun(@(x) x.cachedOrDefault('isUnitary', false), factors);
            factorsAllNonUnitary = cellfun(@(x) ~x.cachedOrDefault('isUnitary', true), factors);
            d = prod(cellfun(@(f) f.dimension, factors));
            args = cell(1, 0);
            if all(factorsAllUnitary)
                args = {'isUnitary' true};
            elseif all(factorsAllNonUnitary)
                args = {'isUnitary' false};
            end
            self@replab.Rep(group, field, d, args{:});
            self.factors = factors;
        end

        function n = nFactors(self)
        % Returns the number of factors in the tensor product
        %
        % Returns:
        %   integer: Number of factors
            n = length(self.factors);
        end

        function f = factor(self, i)
        % Returns a factor in the tensor product
        %
        % Args:
        %   i (integer): Index of a factor
        %
        % Returns:
        %   `+replab.Rep`: Representation corresponding to the i-th factor
            f = self.factors{i};
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rho = image_exact(self, g)
            if self.dimension == 0
                rho = replab.cyclotomic.zeros(0, 0);
            elseif isempty(self.factors)
                rho = replab.cyclotomic.eye(1);
            else
                rho = self.factors{1}.image(g, 'exact');
                for i = 2:self.nFactors
                    rho = kron(rho, self.factors{i}.image(g, 'exact'));
                end
            end
        end

        function rho = image_double_sparse(self, g)
            if self.dimension == 0
                rho = [];
            elseif isempty(self.factors)
                rho = 1;
            else
                rho = self.factors{1}.image(g, 'double/sparse');
                for i = 2:self.nFactors
                    rho = kron(rho, self.factors{i}.image(g, 'double/sparse'));
                end
            end
        end

        function e = computeErrorBound(self)
        % we compute
        % e = [(1 + e1/c1)(1 + e2/c2)...(1 + en/cn) - 1]*(c1 c2 ... cn)
            E = cellfun(@(rep) rep.errorBound, self.factors);
            C = cellfun(@(rep) rep.conditionNumberEstimate, self.factors);
            e = (prod(1 + E./C) - 1)*prod(C);
        end

        function c = computeConditionNumberEstimate(self)
            c = prod(cellfun(@(rep) rep.conditionNumberEstimate, self.factors));
        end

        function rep = computeUnitarize(self)
            srs = cellfun(@(rep) rep.unitarize, self.factors, 'uniform', 0);
            if all(cellfun(@(rep) rep.isExact, srs))
                As = cellfun(@(sr) sr.A('exact'), srs, 'uniform', 0);
                Ainvs = cellfun(@(sr) sr.Ainv('exact'), srs, 'uniform', 0);
                A = replab.numerical.kron(As, 'exact');
                Ainv = replab.numerical.kron(Ainvs, 'exact');
            else
                As = cellfun(@(sr) sr.A('double/sparse'), srs, 'uniform', 0);
                Ainvs = cellfun(@(sr) sr.Ainv('double/sparse'), srs, 'uniform', 0);
                A = replab.numerical.kron(As, 'double/sparse');
                Ainv = replab.numerical.kron(Ainvs, 'double/sparse');
            end
            rep = replab.SimilarRep(self, A, AInv);
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'factors';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:self.nFactors
                names{1, end+1} = sprintf('factor(%d)', i);
                values{1, end+1} = self.factor(i);
            end
        end

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        % Rep

        function b = isExact(self)
            b = all(cellfun(@(f) f.isExact, self.factors));
        end

    end

end
