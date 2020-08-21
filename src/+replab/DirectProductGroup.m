classdef DirectProductGroup < replab.Group
% Describes an external direct product of groups
%
% This is an abstract base class. Use `.DirectProductGroup.make` or `.CompactGroup.directProduct` to construct an instance.

    properties (SetAccess = protected)
        factors % (cell(1,\*) of subclass of `.Group`): Factor groups
    end

    methods (Static)

        function prd = make(factors)
            isFinite = all(cellfun(@(g) isa(g, 'replab.FiniteGroup'), factors));
            if isFinite
                prd = replab.prods.DirectProductOfFiniteGroups(factors);
            else
                prd = replab.prods.DirectProductOfCompactGroups(factors);
            end
        end

    end

    methods

        function self = DirectProductGroup(factors)
        % Constructs a direct product of groups
        %
        % Args:
        %   factors (cell(1,\*) of `.Group`): Factor groups
            self.factors = factors;
            self.identity = cellfun(@(f) f.identity, factors, 'uniform', 0);
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

        function rep = directSumFactorRep(self, factorReps)
        % Constructs a direct sum representation
        %
        % Args:
        %   factorReps (row cell array): Representations for each of the factor groups
        %                                i.e. factorReps{i} is a representation of factor(i)
        %
        % Returns:
        %   `+replab.Rep`: A direct sum representation
            rep = replab.directproduct.SumRep(self, factorReps);
        end

        function rep = directSumFactorRepFun(self, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Args:
        %   fun (function_handle): A function valid for each factor group that maps the group and its index
        %                          to one of its representations
        %
        % Returns:
        %   replab.Rep: A direct sum representation
            reps = arrayfun(@(i) fun(self.factor(i), i), 1:self.nFactors, 'uniform', 0);
            rep = self.directSumFactorRep(reps);
        end

        function rep = tensorFactorRep(self, factorReps)
        % Constructs a tensor product representation
        %
        % Args:
        %   factorReps (row cell array): Representations for each of the factor groups
        %                                i.e. factorReps{i} is a representation of factor(i)
        %
        % Returns:
        %   replab.Rep: A tensor representation
            rep = replab.directproduct.TensorRep(self, factorReps);
        end

        function rep = tensorFactorRepFun(self, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Args:
        %   fun (function_handle): A function valid for each factor group that maps the group and its index
        %                          to one of its representations
        %
        % Returns:
        %   replab.Rep: A tensor representation
            reps = arrayfun(@(i) fun(self.factor(i), i), 1:self.nFactors, 'uniform', 0);
            rep = self.tensorFactorRep(reps);
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
            m = self.factor(i).morphismByFunction(self, @(g) replab.DirectProductGroup.updateCellArray(self.identity, i, g));
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
            m = self.morphismByFunction(self.factor(i), @(g) g{i});
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

    end

end
