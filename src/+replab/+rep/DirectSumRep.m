classdef DirectSumRep < replab.Rep
% A direct sum of representations, such that images are diagonal by blocks

    properties
        factors % (cell(1,\*) of `+replab.Rep`): Contained subrepresentations
    end

    methods

        function self = DirectSumRep(group, field, factors)
        % Constructs a direct sum from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   group (`+replab.CompactGroup`): Common group
        %   field ({'R', 'C'}): Real or complex field
        %   blocks (cell(1,\*) of `+replab.Rep`): Subrepresentations
            replab.rep.assertCompatibleFactors(group, field, factors);
            factorsAllUnitary = cellfun(@(x) x.cachedOrDefault('isUnitary', false), factors);
            factorsAllNonUnitary = cellfun(@(x) ~x.cachedOrDefault('isUnitary', true), factors);
            d = sum(cellfun(@(f) f.dimension, factors));
            args = cell(1, 0);
            if all(factorsAllUnitary)
                args = {'isUnitary', true};
            elseif all(factorsAllNonUnitary)
                args = {'isUnitary', false};
            end
            self@replab.Rep(group, field, d, args{:});
            self.factors = factors;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_someFactorsArePermutationSimilarReps(self, options)
        % Rewrite rule: move permutation similarity transforms before performing the direct sum
            isPermSimilar = cellfun(@(f) isa(f, 'replab.SimilarRep') && f.isPermutation, self.factors);
            if any(isPermSimilar)
                A = sparse(0, 0);
                Ainv = sparse(0, 0);
                n = self.nFactors;
                newFactors = cell(1, n);
                for i = 1:n
                    f = self.factor(i);
                    if isPermSimilar(i)
                        A = blkdiag(A, f.A_internal);
                        Ainv = blkdiag(Ainv, f.Ainv_internal);
                        newFactors{i} = f.parent;
                    else
                        A = blkdiag(A, speye(f.dimension));
                        Ainv = blkdiag(Ainv, speye(f.dimension));
                        newFactors{i} = f;
                    end
                end
                res = replab.SimilarRep(replab.rep.DirectSumRep(self.group, self.field, newFactors), A, Ainv);
            else
                res = [];
            end
        end

        function res = rewriteTerm_removeTrivialFactors(self, options)
        % Rewrite rule: remove trivial factors
            mask = cellfun(@(f) f.dimension == 0, self.factors);
            if any(mask)
                res = replab.rep.DirectSumRep(self.group, self.field, self.factors(~mask));
            else
                res = [];
            end
        end

        function res = rewriteTerm_factorIsDirectSum(self, options)
        % Rewrite rule: if any of the factors is a direct sum itself, collapse the sums
            if any(cellfun(@(f) isa(f, 'replab.rep.DirectSumRep'), self.factors))
                newFactors = cell(1, 0);
                for i = 1:length(self.factors)
                    f = self.factor(i);
                    if isa(f, 'replab.rep.DirectSumRep')
                        newFactors = horzcat(newFactors, f.factors);
                    else
                        newFactors{1,end+1} = f;
                    end
                end
                res = replab.rep.DirectSumRep(self.group, self.field, newFactors);
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

    end

    methods

        function n = nFactors(self)
        % Returns the number of factors in the direct sum
        %
        % Returns:
        %   integer: Number of subrepresentations composing this representation
            n = length(self.factors);
        end

        function f = factor(self, i)
        % Returns a factor in the direct sum
        %
        % Args:
        %   i (integer): Index of a factor/block
        %
        % Returns:
        %   `+replab.Rep`: Representation corresponding to the i-th factor
            f = self.factors{i};
        end

    end

    methods (Access = protected) % Implementations

        function c = decomposeTerm(self)
            c = self.factors;
        end

        function r = composeTerm(self, newFactors)
            r = replab.rep.DirectSumRep(self.group, self.field, newFactors);
        end

        function rho = image_exact(self, g)
            if self.dimension == 0
                rho = replab.cyclotomic.zeros(0, 0);
            else
                rhos = cellfun(@(rep) rep.image(g, 'exact'), self.factors, 'uniform', 0);
                rho = blkdiag(rhos{:});
            end
        end

        function rho = image_double_sparse(self, g)
            rhos = cellfun(@(rep) rep.image(g, 'double/sparse'), self.factors, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end

        function e = computeErrorBound(self)
            e = sum(cellfun(@(rep) rep.errorBound, self.factors));
        end

        function c = computeConditionNumberEstimate(self)
            c = max(cellfun(@(rep) rep.conditionNumberEstimate, self.factors));
        end

        function K = computeKernel(self)
            if self.nFactors == 0
                K = self.group;
            else
                K = self.factor(1).kernel;
                for i = 2:self.nFactors;
                    K = K.intersection(self.factor(i).kernel);
                end
            end
        end

        function rep = computeUnitarize(self)
            srs = cellfun(@(rep) rep.unitarize, self.factors, 'uniform', 0);
            if self.nFactors == 0
                rep = replab.SimilarRep(self, replab.cyclotomic.zeros(0, 0), replab.cyclotomic.zeros(0, 0));
                return
            end
            if all(cellfun(@(rep) rep.isExact, srs))
                As = cellfun(@(sr) sr.A('exact'), srs, 'uniform', 0);
                Ainvs = cellfun(@(sr) sr.Ainv('exact'), srs, 'uniform', 0);
                Ainv = replab.numerical.kron(Ainvs, 'exact');
            else
                As = cellfun(@(sr) sr.A('double/sparse'), srs, 'uniform', 0);
                Ainvs = cellfun(@(sr) sr.Ainv('double/sparse'), srs, 'uniform', 0);
            end
            A = blkdiag(As{:});
            Ainv = blkdiag(Ainvs{:});
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
