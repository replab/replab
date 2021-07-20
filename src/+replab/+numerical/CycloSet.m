classdef CycloSet < replab.Str
% Stores a list of cyclotomic matrices

    properties (SetAccess = protected)
        nRows % (integer): Number of matrix rows
        nCols % (integer): Number of matrix columns
        cycloSet % (cyclo.CycloSet): JVM object storing the matrices
    end

    methods

        function self = CycloSet(nRows, nCols)
            self.nRows = nRows;
            self.nCols = nCols;
            self.cycloSet = javaObject('cyclo.CycloSet', nRows, nCols);
        end

        function s = nElements(self)
            s = self.cycloSet.size;
        end

        function c = at(self, i)
        % Returns the matrix for the given index
            c = replab.cyclotomic(self.cycloSet.at(i-1), [self.nRows self.nCols]);
        end

        function ind = insert(self, c)
        % Inserts the given matrix
            n = self.nElements;
            ind = self.cycloSet.findOrAdd(c.data);
            assert(ind == n, 'Element was already present');
            ind = ind + 1;
        end

        function ind = find(self, c)
            ind = self.cycloSet.find(c.data) + 1;
        end

        function ind = update(self, c)
            ind = self.cycloSet.findOrAdd(c.data);
            ind = ind + 1;
        end

    end

end
