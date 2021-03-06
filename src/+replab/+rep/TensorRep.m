classdef TensorRep < replab.Rep
% A tensor product of representations
%
% All factor representations must be defined on the same group

    properties (SetAccess = protected)
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
        %   factors (cell(1,\*) of `+replab.Rep`): Factor representations
            replab.rep.assertCompatibleFactors(group, field, factors);
            d = prod(cellfun(@(f) f.dimension, factors));
            factorsAreUnitary = cellfun(@(x) x.isUnitary, factors);
            self@replab.Rep(group, field, d, 'isUnitary', all(factorsAreUnitary));
            self.factors = factors;
        end

    end

    methods % Simplification rules

        % TODO
% $$$         function res = rewriteTerm_someFactorsArePermutationSimilarReps(self, options)
% $$$             isPermSimilar = cellfun(@(f) isa(f, 'replab.SimilarRep') && f.isPermutation, self.factors);
% $$$             if any(isPermSimilar)
% $$$                 A = speye(1);
% $$$                 Ainv = speye(1);
% $$$                 n = self.nFactors;
% $$$                 newFactors = cell(1, n);
% $$$                 for i = 1:n
% $$$                     f = self.factor(i);
% $$$                     if isPermSimilar(i)
% $$$                         A = kron(A, f.A_internal);
% $$$                         Ainv = kron(Ainv, f.Ainv_internal);
% $$$                         newFactors{i} = f.parent;
% $$$                     else
% $$$                         A = kron(A, speye(f.dimension));
% $$$                         Ainv = kron(Ainv, speye(f.dimension));
% $$$                         newFactors{i} = f;
% $$$                     end
% $$$                 end
% $$$                 res = replab.SimilarRep(replab.rep.TensorRep(self.group, self.field, newFactors), A, Ainv);
% $$$             else
% $$$                 res = [];
% $$$             end
% $$$         end

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

        function M = matrixRowAction_double_sparse(self, g, M)
            e = size(M, 2);
            n = self.nFactors;
            D = self.dimension;
            d = cellfun(@(f) f.dimension, self.factors);
            % in the Matlab convention, we would have a tensor
            % dn dn-1 ... d1 e
            % first step, we multiply (upper case when has been *)
            % Dn dn-1 ... d1 e
            % then we transpose
            % dn-1 ... d1 e Dn
            % Dn-1 Dn-2 ... d1 e Dn
            % we transpose
            % Dn-2 ... d1 e Dn Dn-1
            % ...
            % d1 e Dn ... D2
            % D1 e Dn ... D2
            % e Dn ... D1
            % last transpose

            % first transpose across Dn | dn-1 ... e
            % second transpose across
            for i = n:-1:1
                F = self.factor(i);
                M = reshape(M, [d(i) e*D/d(i)]);
                M = F.matrixRowAction(g, M);
                M = M.';
            end
            M = reshape(M, [e D]).';
        end

        function M = matrixColAction_double_sparse(self, g, M)
            e = size(M, 1);
            n = self.nFactors;
            D = self.dimension;
            dims = cellfun(@(f) f.dimension, self.factors);
            % in the Matlab convention, we would have a tensor
            % e dn ... d2 d1
            % first step, we multiply (upper case when has been *)
            % e dn ... d2 D1
            % then we transpose
            % D1 e dn ... d2
            %
            % Dn .. D1 e
            for i = 1:n
                F = self.factor(i);
                M = reshape(M, [e*D/dims(i) dims(i)]);
                M = F.matrixColAction(g, M);
                M = M.';
            end
            M = reshape(M, [D e]).';
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
            rep = self.subRep(Ainv, 'projection', A, 'isUnitary', true);
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

        function p = invariantBlocks(self)
            if self.nFactors == 0
                p = replab.Partition.fromBlocks({[1]});
            else
                dimL = self.factor(1).dimension;
                blocksL = self.factor(1).invariantBlocks.blocks;
                for i = 2:self.nFactors
                    dimR = self.factor(i).dimension;
                    blocksR = self.factor(i).invariantBlocks.blocks;
                    res = cell(1, 0);
                    for l = 1:length(blocksL)
                        blkL = blocksL{l};
                        maskL = zeros(1, dimL);
                        maskL(blkL) = 1;
                        for r = 1:length(blocksR)
                            blkR = blocksR{r};
                            maskR = zeros(1, dimR);
                            maskR(blkR) = 1;
                            res{1,end+1} = find(kron(maskL, maskR));
                        end
                    end
                    blocksL = res;
                end
                p = replab.Partition.fromBlocks(blocksL);
            end
        end

        function b = hasTorusImage(self)
            b = all(cellfun(@(f) f.hasTorusImage, self.factors));
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            r = self.group.reconstruction.source.n; % torus rank
            torusMap = zeros(1, r);
            torusInjection = 1;
            torusProjection = 1;
            d = 1;
            for i = 1:self.nFactors
                f = self.factor(i);
                [tm, ti, tp] = f.torusImage;
                torusInjection = kron(torusInjection, ti);
                torusProjection = kron(torusProjection, tp);
                di = f.dimension;
                torusMap1 = zeros(d * di, r);
                for j = 1:r
                    col = bsxfun(@plus, torusMap(:,j)', tm(:,j)); % respects the kron convention
                    torusMap1(:, j) = col(:);
                end
                torusMap = torusMap1;
                d = d * di;
            end
        end

    end

end
