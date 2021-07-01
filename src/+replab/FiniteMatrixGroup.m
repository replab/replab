<<<<<<< HEAD
classdef FiniteMatrixGroup < replab.NiceFiniteGroup

    properties (SetAccess = protected)
        matrixSize % (integer): Matrix size
        elementSet % (`+replab.+numerical.CycloSet`): List of elements (Java object)
    end

    properties (Access = protected)
        permutationImages % (integer(\*,\*)): Each row gives the permutation representation of an element in the regular representation
=======
classdef FiniteMatrixGroup < replab.FiniteGroup

    properties (SetAccess = protected)
        matrixSize % (integer): Matrix size
        elements % (cyclo.CycloSet): List of elements (Java object)
>>>>>>> Starting finite matrix group support
    end

    methods

<<<<<<< HEAD
        function self = FiniteMatrixGroup(matrixSize, generators, varargin)
            generators = cellfun(@(g) replab.cyclotomic(g), generators, 'uniform', 0);
            args = struct('relators', []);
            identity = replab.cyclotomic.eye(matrixSize);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            type = 'self';
            self@replab.NiceFiniteGroup(identity, generators, type, restArgs{:});
            self.matrixSize = matrixSize;
            self.elementSet = replab.matrix.dimino(generators, identity);
            order = self.elementSet.nElements;
            self.permutationImages = zeros(order, order);
=======
        function self = FiniteMatrixGroup(matrixSize, generators)
            self.identity = replab.cyclotomic.eye(matrixSize);
>>>>>>> Starting finite matrix group support
        end

    end

    methods % Implementations

<<<<<<< HEAD
        % Domain

        function l = eqv(self, x, y)
            l = full(all(all(x == y)));
        end

        function g = sample(self)
            g = self.elementSet.at(randi(self.elementSet.nElements));
        end

        % Monoid

        function z = compose(self, x, y)
            z = x * y;
        end

        % Group

=======
>>>>>>> Starting finite matrix group support
        function xInv = inverse(self, x)
            xInv = inv(x);
        end

        % FiniteGroup

        function G = withGeneratorNames(self, newNames)
            if isequal(self.generatorNames, newNames)
                G = self;
                return
            end
            G = replab.FiniteMatrixGroup(self.matrixSize, self.generators, 'generatorNames', newNames, 'type', self.type);
        end

        % NiceFiniteGroup

        function res = hasSameTypeAs(self, rhs)
            res = self.type.id == rhs.type.id;
        end

        function p = niceImage(self, g)
            i = self.elementSet.find(g);
            assert(i > 0);
            if all(self.permutationImages(i, :) == 0)
                for j = 1:self.elementSet.nElements
                    self.permutationImages(i, j) = self.elementSet.find(g * self.elementSet.at(j));
                end
            end
            p = self.permutationImages(i, :);
        end

    end

end
