classdef DirectProductGroup < replab.CompactGroup
% Describes an external direct product of compact groups
%
% This is an abstract base class. Use `.DirectProductGroup.make` or `.CompactGroup.directProduct` to construct an instance.
%
% Constructors are defined in subclasses, as for other group products.

    properties (SetAccess = protected)
        factors % (cell(1,\*) of subclass of `.CompactGroup`): Factor groups
    end

    methods (Static)

        function prd = make(factors)
            isCompact = cellfun(@(g) isa(g, 'replab.CompactGroup'), factors);
            assert(all(isCompact), 'All factors must be compact');
            isFinite = cellfun(@(g) isa(g, 'replab.FiniteGroup'), factors);
            if all(isFinite)
                prd = replab.prods.DirectProductGroup_finite(factors);
            else
                prd = replab.prods.DirectProductGroup_compact(factors);
            end
        end

    end

    methods % Factor manipulation

        function n = nFactors(self)
        % Returns the number of factors in this direct product
        %
        % Returns:
        %   integer: Number of factors
            n = length(self.factors);
        end

        function f = factor(self, i)
        % Returns a factor in this direct product
        %
        % Args:
        %   i (integer): Factor idnex
        %
        % Returns:
        %   `.Group`: The ``i``-th factor
            f = self.factors{i};
        end

    end

    methods % Representations

        function rep = directSumFactorRep(self, field, factorReps)
        % Constructs a direct sum representation
        %
        % Args:
        %   field ({'R', 'C'}): Field
        %   factorReps (row cell array): Representations for each of the factor groups
        %                                i.e. factorReps{i} is a representation of factor(i)
        %
        % Returns:
        %   `+replab.Rep`: A direct sum representation
            factors = arrayfun(@(i) self.projection(i).andThen(factorReps{i}), 1:self.nFactors, 'uniform', 0);
            rep = self.directSumRep(field, factors);
        end

        function rep = directSumFactorRepFun(self, field, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Args:
        %   field ({'R', 'C'}): Field
        %   fun (function_handle): A function valid for each factor group that maps the group and its index
        %                          to one of its representations
        %
        % Returns:
        %   `+replab.Rep`: A direct sum representation
            reps = arrayfun(@(i) fun(self.factor(i), i), 1:self.nFactors, 'uniform', 0);
            rep = self.directSumFactorRep(field, reps);
        end

        function rep = tensorFactorRep(self, field, factorReps)
        % Constructs a tensor product representation
        %
        % Args:
        %   field ({'R', 'C'}): Field
        %   factorReps (row cell array): Representations for each of the factor groups
        %                                i.e. factorReps{i} is a representation of factor(i)
        %
        % Returns:
        %   `+replab.Rep`: A tensor representation
            factors = arrayfun(@(i) self.projection(i).andThen(factorReps{i}), 1:self.nFactors, 'uniform', 0);
            rep = self.tensorRep(field, factors);
        end

        function rep = tensorFactorRepFun(self, field, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Args:
        %   field ({'R', 'C'}): Field
        %   fun (function_handle): A function valid for each factor group that maps the group and its index
        %                          to one of its representations
        %
        % Returns:
        %   replab.Rep: A tensor representation
            reps = arrayfun(@(i) fun(self.factor(i), i), 1:self.nFactors, 'uniform', 0);
            rep = self.tensorFactorRep(field, reps);
        end

    end

    methods % Morphisms

        function m = embedding(self, i)
        % Returns the morphism embedding the i-th factor into the direct product
        %
        % Example:
        %   >>> S2 = replab.S(2);
        %   >>> D = S2.directProduct(S2);
        %   >>> m = D.embedding(1);
        %   >>> D.eqv({[2 1] [1 2]}, m.imageElement([2 1]))
        %       1
        %
        % Args:
        %   i (integer): Factor index
        %
        % Returns:
        %   `.Morphism`: The embedding
            p = self.torusBlocks;
            b = p.block(i);
            tm = full(sparse(1:length(b), b, ones(1, length(b)), length(b), p.n));
            m = self.factor(i).morphismByFunction(self, @(g) replab.DirectProductGroup.updateCellArray(self.identity, i, g), tm);
        end

        function m = projection(self, i)
        % Returns the morphism projecting this group into its i-th factor
        %
        % Example:
        %   >>> S2 = replab.S(2);
        %   >>> D = S2.directProduct(S2);
        %   >>> m = D.projection(1);
        %   >>> S2.eqv([2 1], m.imageElement({[2 1] [1 2]}))
        %       1
        %
        % Args:
        %   i (integer): Factor index
        %
        % Returns:
        %   `.Morphism`: The projection
            p = self.torusBlocks;
            b = p.block(i);
            tm = full(sparse(b, 1:length(b), ones(1, length(b)), p.n, length(b)));
            m = self.morphismByFunction(self.factor(i), @(g) g{i}, tm);
        end

    end

    methods (Static, Access = protected)

        function c = updateCellArray(c, i, v)
        % Updates the i-th element of a cell array
            c{i} = v;
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Group(self);
            names{1, end+1} = 'factors';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nFactors
                names{1, end+1} = sprintf('factor(%d)', i);
                values{1, end+1} = self.factor(i);
            end
        end

        function s = headerStr(self)
            if self.inCache('order')
                s = sprintf('Direct product group with %d factors of order %s', self.nFactors, strtrim(num2str(self.order)));
            else
                s = sprintf('Direct product group with %d factors', self.nFactors);
            end
        end

        % Domain

        function b = eqv(self, x, y)
            b = true;
            for i = 1:self.nFactors
                if ~self.factor(i).eqv(x{i}, y{i})
                    b = false;
                    return
                end
            end
        end

        function g = sample(self)
            g = cell(1, self.nFactors);
            for i = 1:self.nFactors
                g{i} = self.factor(i).sample;
            end
        end

        % Monoid

        function z = compose(self, x, y)
            z = cell(1, self.nFactors);
            for i = 1:self.nFactors
                z{i} = self.factor(i).compose(x{i}, y{i});
            end
        end

        % Group

        function xInv = inverse(self, x)
            xInv = cell(1, self.nFactors);
            for i = 1:self.nFactors
                xInv{i} = self.factor(i).inverse(x{i});
            end
        end

        % Compact group

        function p = torusBlocks(self)
            b = all(cellfun(@(f) f.hasReconstruction, self.factors));
            if ~b
                p = [];
                return
            end
            n = self.nFactors;
            dims = cellfun(@(f) f.reconstruction.source.n, self.factors);
            blocks = cell(1, n);
            shift = 0;
            for i = 1:n
                blocks{1,i} = shift+(1:dims(i));
                shift = shift + dims(i);
            end
            p = replab.Partition.fromBlocks(blocks);
        end

        function b = hasReconstruction(self)
            b = ~isempty(self.torusBlocks);
        end

        function [mu, R] = reconstruction(self)
            p = self.torusBlocks;
            n = self.nFactors;
            f = @(t) arrayfun(@(i) self.factor(i).reconstruction.imageElement(t(p.block(i))), 1:n, 'uniform', 0);
            T = replab.TorusGroup(p.n);
            sets = cell(1, 0);
            for i = 1:n
                [~, Ri] = self.factor(i).reconstruction;
                e = self.embedding(i);
                sets = horzcat(sets, cellfun(@(S) cellfun(@(s) e.imageElement(s), S, 'uniform', 0), Ri.sets, 'uniform', 0));
            end
            R = replab.SetProduct(self, sets, true);
            mu = T.morphismByFunction(self, f);
        end

    end

end
