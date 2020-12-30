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
            d = prod(cellfun(@(f) f.dimension, factors));
            factorsAllUnitary = cellfun(@(x) x.knownUnitary, factors);
            factorsAllNonUnitary = cellfun(@(x) x.knownNonUnitary, factors);
            if all(factorsAllUnitary)
                args = {'isUnitary' true};
            elseif all(factorsAllNonUnitary)
                args = {'isUnitary' false};
            else
                args = cell(1, 0);
            end
            self@replab.Rep(group, field, d, args{:});
            self.factors = factors;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_someFactorsArePermutationSimilarReps(self, options)
            isPermSimilar = cellfun(@(f) isa(f, 'replab.SimilarRep') && f.isPermutation, self.factors);
            if any(isPermSimilar)
                A = speye(1);
                Ainv = speye(1);
                n = self.nFactors;
                newFactors = cell(1, n);
                for i = 1:n
                    f = self.factor(i);
                    if isPermSimilar(i)
                        A = kron(A, f.A_internal);
                        Ainv = kron(Ainv, f.Ainv_internal);
                        newFactors{i} = f.parent;
                    else
                        A = kron(A, speye(f.dimension));
                        Ainv = kron(Ainv, speye(f.dimension));
                        newFactors{i} = f;
                    end
                end
                res = replab.SimilarRep(replab.rep.TensorRep(self.group, self.field, newFactors), A, Ainv);
            else
                res = [];
            end
        end

        function res = rewriteTerm_moveDirectSumToFirstFactor(self, options)
        % If a factor is a direct sum, move it as the first factor
            i = find(cellfun(@(f) isa(f, 'replab.rep.DirectSumRep'), self.factors), 1);
            if ~isempty(i)
                n = self.nFactors;
                perm = 1:n;
                perm([(n-i+1) n]) = [n (n-i+1)];
                D = self.dimension;
                dims = cellfun(@(f) f.dimension, self.factors);
                ind = reshape(1:D, fliplr(dims));
                ind = permute(ind, perm);
                ind = ind(:)';
                Ainv = sparse(1:D, ind, ones(1, D), D, D);
                A = sparse(ind, 1:D, ones(1, D), D, D);
                newFactors = self.factors;
                newFactors{1} = self.factor(i);
                newFactors{i} = self.factor(1);
                newTP = replab.rep.TensorRep(self.group, self.field, newFactors);
                res = replab.SimilarRep(newTP, A, Ainv);
            else
                res = [];
            end
        end

        function res = rewriteTerm_firstFactorIsDirectSum(self, options)
        % If the first factor is a direct sum, distribute the tensor product
            if self.nFactors > 0 && isa(self.factor(1), 'replab.rep.DirectSumRep')
                f1 = self.factor(1);
                n = f1.nFactors;
                terms = cell(1, n);
                for i = 1:n
                    newFactors = self.factors;
                    newFactors{1} = f1.factor(i);
                    terms{i} = replab.rep.TensorRep(self.group, self.field, newFactors);
                end
                res = replab.rep.DirectSumRep(self.group, self.field, terms);
            else
                res = [];
            end
        end

        function res = rewriteTerm_factorIsTrivial(self, options)
            mask = cellfun(@(f) f.dimension == 1 && f.cachedOrDefault('trivialDimension', -1) == 1, self.factors);
            if any(mask)
                res = replab.rep.TensorRep(self.group, self.field, self.factors(~mask));
            else
                res = [];
            end
        end

        function res = rewriteTerm_hasOneFactor(self, options)
        % Rewrite rule: removes the direct sum if it has a single factor
            if self.nFactors == 1
                res = self.factor(1);
            else
                res = [];
            end
        end

        function res = rewriteTerm_factorIsTensor(self, options)
        % If any of the factors is a tensor product itself, collapse the products
            if any(cellfun(@(f) isa(f, 'replab.rep.TensorRep'), self.factors))
                newFactors = cell(1, 0);
                for i = 1:length(self.factors)
                    f = self.factor(i);
                    if isa(f, 'replab.rep.TensorRep')
                        newFactors = horzcat(newFactors, f.factors);
                    else
                        newFactors{1,end+1} = f;
                    end
                end
                res = replab.rep.TensorRep(self.group, self.field, newFactors);
            else
                res = [];
            end
        end

    end

    methods

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

        function rep = computeDouble(self)
            rep = replab.rep.TensorRep(self.group, self.field, cellfun(@(f) double(f), self.factors));
        end

        function c = decomposeTerm(self)
            c = self.factors;
        end

        function r = composeTerm(self, newFactors)
            r = replab.rep.TensorRep(self.group, self.field, newFactors);
        end

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
            E = cellfun(@(rep) rep.errorBound, self.factors);
            C = cellfun(@(rep) rep.conditionNumberEstimate, self.factors);
            % we compute
            % e = e1*c2*...*cn + c1*e2*...*cn + ...
            e = 0;
            for i = 1:self.nFactors
                e = e + E(i)*prod(C(1:i-1))*prod(C(i+1:end));
            end
        end

        function c = computeConditionNumberEstimate(self)
            c = prod(cellfun(@(rep) rep.conditionNumberEstimate, self.factors));
        end

        function b = computeIsUnitary(self)
            b = all(cellfun(@(r) r.isUnitary, self.factors));
        end

        function rep = computeUnitarize(self)
            srs = cellfun(@(rep) rep.unitarize, self.factors, 'uniform', 0);
            if all(cellfun(@(rep) rep.isExact, srs))
                As = cellfun(@(sr) sr.A('exact'), srs, 'uniform', 0);
                Ainvs = cellfun(@(sr) sr.Ainv('exact'), srs, 'uniform', 0);
                A = replab.numerical.multikron(As, 'exact');
                Ainv = replab.numerical.multikron(Ainvs, 'exact');
            else
                As = cellfun(@(sr) sr.A('double/sparse'), srs, 'uniform', 0);
                Ainvs = cellfun(@(sr) sr.Ainv('double/sparse'), srs, 'uniform', 0);
                A = replab.numerical.multikron(As, 'double/sparse');
                Ainv = replab.numerical.multikron(Ainvs, 'double/sparse');
            end
            rep = replab.SimilarRep(self, A, Ainv);
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
